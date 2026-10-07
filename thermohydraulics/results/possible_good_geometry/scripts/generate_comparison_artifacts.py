#!/usr/bin/env python3
"""
generate_comparison_artifacts.py

Generates publication-quality comparative figures (300 DPI) and tables:
Comparing GEOM_OPT15_EXTENDED vs SIMULATION 4B REFERENCE.
Highlights:
- Operational lifespan extension (>15 years vs 12.5 years)
- Thermal margins (Clad creep margin, Fuel dissociation margin, Moderator < 950 K)
- Control drum worth & shutdown margins
- Mass budget (< 600 kg) & pressure drop (< 35 kPa)
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14,
    "lines.linewidth": 2.0,
    "lines.markersize": 6,
    "grid.alpha": 0.35,
    "grid.linestyle": "--"
})

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/possible_good_geometry"
DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Load depletion datasets
dep_opt15_path = os.path.join(DATA_DIR, "opt15_depletion_2.5yr_steps.csv")
dep_4b_path = "/home/aleri/OpenMC_projects/first_test/final_calculation/simulation_fourth/sim4b_al2o3/data/depletion_2.5yr_steps_4b.csv"

def make_depletion_comparison_plot():
    if not os.path.exists(dep_opt15_path) or not os.path.exists(dep_4b_path):
        print("Depletion data not found.")
        return

    df_15 = pd.read_csv(dep_opt15_path)
    df_4b = pd.read_csv(dep_4b_path)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Panel 1: Multiplication factor keff(0 deg) vs Time
    ax1.plot(df_4b["time_years"], df_4b["k_eff_drums_out_0deg"], 's-', color="#d95f02", label="Sim 4B (H=42.45 cm, 124 kgU)")
    ax1.plot(df_15["time_years"], df_15["k_eff_drums_out_0deg"], 'o-', color="#1b9e77", label=r"GEOM_OPT15 (H=45.50 cm, 135 kgU)")
    ax1.axhline(1.000, color="black", linestyle="--", linewidth=1.5, label=r"Critical Threshold ($k_{\mathrm{eff}}=1.000$)")
    ax1.axvline(15.0, color="#7570b3", linestyle=":", linewidth=2.0, label="Target Lifespan (15.0 yr)")
    
    ax1.fill_between([15.0, 30.0], [1.0, 1.0], [1.03, 1.03], color="#1b9e77", alpha=0.08, label=r"Target Zone ($t \geq 15$ yr, $k \geq 1$)")

    ax1.set_xlabel("Operational Time [Years]")
    ax1.set_ylabel(r"Multiplication Factor $k_{\mathrm{eff}}$ ($\theta = 0^\circ$)")
    ax1.set_title(r"(a) Core Reactivity Trajectory & Lifetime")
    ax1.grid(True)
    ax1.legend(loc="upper right")
    ax1.set_xlim(0, 28)
    ax1.set_ylim(0.990, 1.025)

    # Panel 2: Critical Drum Angle theta_crit(t)
    ax2.plot(df_4b["time_years"], df_4b["critical_drum_angle_deg"], 's-', color="#d95f02", label="Sim 4B Critical Angle")
    ax2.plot(df_15["time_years"], df_15["critical_drum_angle_deg"], 'o-', color="#1b9e77", label="GEOM_OPT15 Critical Angle")
    ax2.axhline(0.0, color="red", linestyle="--", linewidth=1.5, label=r"Drum Exhaustion ($\theta_{\mathrm{crit}}=0^\circ$)")
    ax2.axvline(15.0, color="#7570b3", linestyle=":", linewidth=2.0, label="Target Lifespan (15.0 yr)")

    ax2.set_xlabel("Operational Time [Years]")
    ax2.set_ylabel(r"Critical Control Drum Angle $\theta_{\mathrm{crit}}$ [degrees]")
    ax2.set_title(r"(b) Active Control Reserve & Drum Position")
    ax2.grid(True)
    ax2.legend(loc="upper right")
    ax2.set_xlim(0, 28)
    ax2.set_ylim(-5, 90)

    plt.tight_layout()
    out_fig = os.path.join(FIG_DIR, "fig1_lifetime_and_drum_trajectory.png")
    plt.savefig(out_fig, dpi=300)
    plt.close()
    print(f"[+] Depletion figure saved to {out_fig}")


def make_multiphysics_summary_plot():
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Bar chart 1: Core Lifespan
    bars = axes[0].bar(["Sim 4B", "GEOM_OPT15"], [12.50, 25.41], color=["#d95f02", "#1b9e77"], width=0.55)
    axes[0].axhline(15.0, color="red", linestyle="--", linewidth=1.8, label="Target Limit: > 15.0 yr")
    axes[0].set_ylabel("Operational Lifespan [Years]")
    axes[0].set_title(r"(a) Reactivity Lifespan ($k_{\mathrm{eff}} \geq 1.0$)")
    axes[0].grid(axis="y")
    axes[0].legend(loc="upper left")
    axes[0].set_ylim(0, 32)
    for bar in bars:
        yval = bar.get_height()
        axes[0].text(bar.get_x() + bar.get_width()/2.0, yval + 0.8, f"{yval:.1f} yr", ha='center', va='bottom', fontweight='bold')

    # Bar chart 2: Total Mass
    bars2 = axes[1].bar(["Sim 4B", "GEOM_OPT15"], [537.8, 550.5], color=["#d95f02", "#1b9e77"], width=0.55)
    axes[1].axhline(600.0, color="red", linestyle="--", linewidth=1.8, label="Constraint Limit: < 600 kg")
    axes[1].set_ylabel("Total Reactor Subassembly Mass [kg]")
    axes[1].set_title("(b) Reactor Mass Budget")
    axes[1].grid(axis="y")
    axes[1].legend(loc="upper left")
    axes[1].set_ylim(400, 650)
    for bar in bars2:
        yval = bar.get_height()
        axes[1].text(bar.get_x() + bar.get_width()/2.0, yval + 5.0, f"{yval:.1f} kg", ha='center', va='bottom', fontweight='bold')

    # Bar chart 3: Core Pressure Drop
    bars3 = axes[2].bar(["Sim 4B", "GEOM_OPT15"], [19.21, 19.45], color=["#d95f02", "#1b9e77"], width=0.55)
    axes[2].axhline(35.0, color="red", linestyle="--", linewidth=1.8, label="Pumping Limit: < 35 kPa")
    axes[2].set_ylabel(r"Core Pressure Drop $\Delta P$ [kPa]")
    axes[2].set_title("(c) Thermal-Hydraulic Pressure Loss")
    axes[2].grid(axis="y")
    axes[2].legend(loc="upper left")
    axes[2].set_ylim(0, 42)
    for bar in bars3:
        yval = bar.get_height()
        axes[2].text(bar.get_x() + bar.get_width()/2.0, yval + 0.8, f"{yval:.2f} kPa", ha='center', va='bottom', fontweight='bold')

    plt.tight_layout()
    out_fig2 = os.path.join(FIG_DIR, "fig2_multiphysics_constraints_comparison.png")
    plt.savefig(out_fig2, dpi=300)
    plt.close()
    print(f"[+] Multiphysics constraints figure saved to {out_fig2}")


def make_thermal_comparison_plot():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.2))

    materials = ["Fuel Centerline", "Cladding Peak (1D)", "Cladding Peak (2D)", "YH2 Moderator"]
    t_4b = [1246.1, 1218.6, 1194.4, 940.0]
    t_15 = [1246.9, 1214.6, 1191.1, 939.5]

    x = np.arange(len(materials))
    w = 0.35

    rects1 = ax1.bar(x - w/2, t_4b, w, label="Sim 4B", color="#d95f02")
    rects2 = ax1.bar(x + w/2, t_15, w, label="GEOM_OPT15", color="#1b9e77")

    ax1.axhline(1350.0, color="red", linestyle="--", label="Clad Creep Limit (1350 K)")
    ax1.axhline(950.0, color="blue", linestyle=":", label="YH2 Moderator Limit (950 K)")

    ax1.set_ylabel("Operating Temperature [K]")
    ax1.set_title("(a) Peak Component Operating Temperatures")
    ax1.set_xticks(x)
    ax1.set_xticklabels(materials, rotation=15)
    ax1.grid(axis="y")
    ax1.legend(loc="upper right")
    ax1.set_ylim(800, 1450)

    # Margins
    margins = ["Clad Creep (1D)", "Clad Creep (2D)", "Fuel Dissociation", "Moderator Margin"]
    m_4b = [1350.0 - 1218.6, 1350.0 - 1194.4, 2073.0 - 1246.1, 950.0 - 940.0]
    m_15 = [1350.0 - 1214.6, 1350.0 - 1191.1, 2073.0 - 1246.9, 950.0 - 939.5]

    xm = np.arange(len(margins))
    ax2.bar(xm - w/2, m_4b, w, label="Sim 4B", color="#d95f02")
    ax2.bar(xm + w/2, m_15, w, label="GEOM_OPT15", color="#1b9e77")
    ax2.set_ylabel("Safety Thermal Margin [K]")
    ax2.set_title("(b) Design Margins to Material Limits")
    ax2.set_xticks(xm)
    ax2.set_xticklabels(margins, rotation=15)
    ax2.grid(axis="y")
    ax2.legend(loc="upper right")
    ax2.set_ylim(0, 950)

    plt.tight_layout()
    out_fig3 = os.path.join(FIG_DIR, "fig3_thermal_margins_comparison.png")
    plt.savefig(out_fig3, dpi=300)
    plt.close()
    print(f"[+] Thermal margins figure saved to {out_fig3}")


if __name__ == "__main__":
    make_multiphysics_summary_plot()
    make_depletion_comparison_plot()
    make_thermal_comparison_plot()
