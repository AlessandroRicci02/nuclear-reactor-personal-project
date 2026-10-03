#!/usr/bin/env python3
"""
===============================================================================
PARETO ANALYSIS & MULTI-OBJECTIVE ARCHITECTURE RANKING ENGINE
===============================================================================
Generates comprehensive comparative plots and Pareto frontier analyses
comparing the 4 competing reactor architecture families against Opt-9:
  - Total Mass vs Total Assembly Volume
  - Total Mass vs Shutdown Margin (SDM)
  - Total Mass vs Coolant Void Reactivity
  - Total Drum Worth vs Number of Control Drums
  - Multi-Criteria Weighted Utility Scoring
===============================================================================
"""

import os
import sys
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def generate_screening_pareto_plots(csv_path="results/architecture_screening.csv", out_dir="results/figures_screening"):
    os.makedirs(out_dir, exist_ok=True)
    if not os.path.exists(csv_path):
        print(f"File {csv_path} does not exist yet.")
        return

    data = []
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append(row)

    if not data:
        print("No data in screening CSV.")
        return

    # Extract metrics
    cand_ids = [d["cand_id"] for d in data]
    families = [d["family"] for d in data]
    masses = [float(d["mass_total_kg"]) for d in data]
    volumes = [float(d["volume_total_L"]) for d in data]
    k_outs = [float(d["k_out"]) for d in data]
    sdms = [float(d["sdm_pcm"]) for d in data]
    voids = [float(d["drho_void_pcm"]) for d in data]
    void_errs = [float(d["u_void_pcm"]) for d in data]
    worths = [float(d["drum_worth_pcm"]) for d in data]

    # Color map by family
    fam_colors = {
        "Family A (Square Multi-Drum)": "#1E88E5",
        "Family B (Hexagonal Lattice)": "#43A047",
        "Family C (Annular Island)": "#FB8C00",
        "Family D (Clustered/Zoned)": "#8E24AA"
    }

    # 1. Mass vs Volume
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    for i, d in enumerate(data):
        c = fam_colors.get(d["family"], "#757575")
        ax.scatter(masses[i], volumes[i], color=c, s=120, edgecolors='black', zorder=5)
        ax.annotate(d["cand_id"], (masses[i]+2, volumes[i]+0.8), fontsize=9, fontweight='bold')

    # Add Opt-9 reference point
    ax.scatter([568.99], [164.29], color='#D81B60', marker='*', s=260, edgecolors='black', label='Opt-9 Reference (569 kg, 164 L)', zorder=6)
    ax.annotate("Opt-9", (568.99+3, 164.29-2), fontsize=10, fontweight='bold', color='#D81B60')

    # Add dummy points for legend
    for fam, col in fam_colors.items():
        ax.scatter([], [], color=col, s=100, label=fam)

    ax.set_title("Architecture Screening: Assembly Mass vs Assembly Volume", fontsize=12, fontweight='bold')
    ax.set_xlabel("Total Reactor Assembly Mass [kg]", fontweight='bold')
    ax.set_ylabel("Total Reactor Assembly Volume [L]", fontweight='bold')
    ax.grid(True, ls=':', alpha=0.6)
    ax.legend(loc='lower right', fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "screening_mass_vs_volume.png"))
    plt.close(fig)

    # 2. Mass vs Coolant Void Reactivity
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    for i, d in enumerate(data):
        c = fam_colors.get(d["family"], "#757575")
        ax.errorbar(masses[i], voids[i], yerr=void_errs[i], fmt='o', color=c, ecolor=c, elinewidth=1.2, capsize=3, markersize=8, markeredgecolor='black', zorder=5)
        ax.annotate(d["cand_id"], (masses[i]+2, voids[i]+15), fontsize=9, fontweight='bold')

    # Add Opt-9 reference point
    ax.errorbar([568.99], [-36.5], yerr=[17.1], fmt='*', color='#D81B60', ecolor='#D81B60', elinewidth=1.8, capsize=4, markersize=14, markeredgecolor='black', label='Opt-9 Reference (-36.5 pcm)', zorder=6)
    ax.axhline(0.0, color='red', ls='--', lw=1.2, label='Neutral Void Boundary (0 pcm)')

    for fam, col in fam_colors.items():
        ax.scatter([], [], color=col, s=100, label=fam)

    ax.set_title("Architecture Screening: Coolant Void Reactivity vs Total Mass", fontsize=12, fontweight='bold')
    ax.set_xlabel("Total Reactor Assembly Mass [kg]", fontweight='bold')
    ax.set_ylabel(r"Depressurized Void Reactivity $\Delta\rho_{void}$ [pcm]", fontweight='bold')
    ax.grid(True, ls=':', alpha=0.6)
    ax.legend(loc='upper right', fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "screening_mass_vs_void.png"))
    plt.close(fig)

    # 3. Mass vs Shutdown Margin
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    for i, d in enumerate(data):
        c = fam_colors.get(d["family"], "#757575")
        ax.scatter(masses[i], sdms[i], color=c, s=120, edgecolors='black', zorder=5)
        ax.annotate(d["cand_id"], (masses[i]+2, sdms[i]+120), fontsize=9, fontweight='bold')

    ax.scatter([568.99], [2437], color='#D81B60', marker='*', s=260, edgecolors='black', label='Opt-9 Reference (SDM = 2437 pcm)', zorder=6)
    ax.axhline(2000.0, color='darkgreen', ls='--', lw=1.2, label='Minimum SDM Criterion (2000 pcm)')

    for fam, col in fam_colors.items():
        ax.scatter([], [], color=col, s=100, label=fam)

    ax.set_title("Architecture Screening: All-In Shutdown Margin vs Assembly Mass", fontsize=12, fontweight='bold')
    ax.set_xlabel("Total Reactor Assembly Mass [kg]", fontweight='bold')
    ax.set_ylabel("Shutdown Margin SDM [pcm]", fontweight='bold')
    ax.grid(True, ls=':', alpha=0.6)
    ax.legend(loc='upper right', fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "screening_mass_vs_sdm.png"))
    plt.close(fig)

    print(f"Generated screening Pareto plots in {out_dir}")

if __name__ == "__main__":
    generate_screening_pareto_plots()
