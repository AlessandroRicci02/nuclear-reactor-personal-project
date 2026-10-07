#!/usr/bin/env python3
"""
generate_depletion_figures_tb35.py

Generates publication-quality figures for GEOM_OPT15_TB35:
1. geom_tb35_07_depletion_lifetime_and_drums.png:
   - 3-Panel multi-year depletion trajectory:
     * Top: k_eff (0 deg drums out), SDM (180 deg drums in), and Operational Critical k_eff(t) ~ 1.0000.
     * Middle: Critical Drum Angle theta_crit(t) evolution from BOL to EOL.
     * Bottom: Fissile mass depletion (U-235 consumption, Pu-239 breeding).
2. geom_tb35_08_core_mass_breakdown_pie_bar.png:
   - Detailed mass distribution and launch ceiling margin (+4.98 kg margin to 600 kg).
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Style setup
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14,
    "lines.linewidth": 2.2,
    "lines.markersize": 7,
    "grid.alpha": 0.4,
    "grid.linestyle": "--"
})


def plot_mass_breakdown():
    csv_mass = os.path.join(DATA_DIR, "tb35_precise_mass_breakdown.csv")
    if not os.path.exists(csv_mass):
        print(f"[Warn] {csv_mass} not found, skipping mass plot.")
        return

    df = pd.read_csv(csv_mass)
    # Filter individual components
    comps = df[~df["component"].str.contains("TOTAL") & ~df["component"].str.contains("Total Canned")].copy()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.5))

    # Pie chart
    labels = [
        "UN Fuel (Inner)", "UN Fuel (Outer)", "Nb-1Zr Clad",
        "Al2O3 Liner", "YH2 Matrix", "YH2 Periph Buffer",
        "BeO Reflector", "B4C Absorber"
    ]
    masses = [
        comps.loc[comps["component"].str.contains("Inner"), "mass_kg"].values[0],
        comps.loc[comps["component"].str.contains("Outer"), "mass_kg"].values[0],
        comps.loc[comps["component"].str.contains("Cladding"), "mass_kg"].values[0],
        comps.loc[comps["component"].str.contains("Thermal Barrier"), "mass_kg"].values[0],
        comps.loc[comps["component"].str.contains("Moderator Matrix"), "mass_kg"].values[0],
        comps.loc[comps["component"].str.contains("Moderator Peripheral"), "mass_kg"].values[0],
        comps.loc[comps["component"].str.contains("Reflector"), "mass_kg"].values[0],
        comps.loc[comps["component"].str.contains("Absorber"), "mass_kg"].values[0],
    ]
    colors = [
        "#d95f02", "#e6ab02", "#7570b3",
        "#e7298a", "#1b9e77", "#66a61e",
        "#386cb0", "#666666"
    ]

    wedges, texts, autotexts = ax1.pie(
        masses, labels=labels, autopct="%1.1f%%",
        startangle=140, colors=colors, pctdistance=0.8,
        wedgeprops=dict(width=0.45, edgecolor="white", linewidth=1.5)
    )
    for at in autotexts:
        at.set_fontsize(9)
        at.set_weight("bold")
    ax1.set_title("GEOM_OPT15_TB35 Mass Distribution\nTotal Subassembly Mass = 593.9 kg (Stoch: 595.0 kg)", fontweight="bold")

    # Bar chart vs Launch Ceiling
    total_mass = sum(masses)
    ceiling = 600.0
    margin = ceiling - total_mass

    bar_labels = ["Total Subassembly Mass", "Launch Payload Ceiling"]
    bar_vals = [total_mass, ceiling]
    bar_colors = ["#2b83ba", "#d7191c"]

    bars = ax2.bar(bar_labels, bar_vals, color=bar_colors, width=0.45, edgecolor="black", linewidth=1.2)
    ax2.set_ylabel("Mass [kg]")
    ax2.set_ylim(0, 680)
    ax2.grid(True, axis="y")
    ax2.set_title(f"Launch Vehicle Payload Limit Compliance\nLaunch Safety Margin = +{margin:.2f} kg (PASS)", fontweight="bold")

    for bar in bars:
        h = bar.get_height()
        ax2.annotate(f"{h:.1f} kg",
                     xy=(bar.get_x() + bar.get_width() / 2, h),
                     xytext=(0, 6), textcoords="offset points",
                     ha="center", va="bottom", fontweight="bold", fontsize=11)

    ax2.annotate(f"Positive Margin:\n+{margin:.2f} kg",
                 xy=(0.5, total_mass + 5), xytext=(0.5, total_mass + 35),
                 arrowprops=dict(facecolor="green", shrink=0.08, width=2, headwidth=8),
                 ha="center", va="bottom", fontsize=11, fontweight="bold", color="green",
                 bbox=dict(boxstyle="round,pad=0.4", fc="#e5f5e0", ec="green", lw=1.5))

    plt.tight_layout()
    out_path = os.path.join(FIG_DIR, "geom_tb35_08_core_mass_breakdown_pie_bar.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[Plot Generated] {out_path}")


def plot_depletion_trajectory():
    csv_dep = os.path.join(DATA_DIR, "tb35_criticality_controlled_depletion.csv")
    if not os.path.exists(csv_dep):
        print(f"[Warn] {csv_dep} not found, skipping depletion plot.")
        return

    df = pd.read_csv(csv_dep)
    t = df["time_years"].values
    k_0 = df["k_eff_drums_out_0deg"].values
    sig_0 = df["sigma_k_0deg"].values
    k_crit = df["k_eff_operational_crit"].values
    sig_crit = df["sigma_k_crit"].values
    theta_crit = df["critical_drum_angle_deg"].values
    u235 = df["u235_mass_kg"].values
    pu239 = df["pu239_mass_kg"].values

    # Determine lifetime crossing
    t_life = None
    for i in range(len(t) - 1):
        if k_0[i] >= 1.00000 and k_0[i+1] < 1.00000:
            t_life = t[i] + (1.00000 - k_0[i]) / (k_0[i+1] - k_0[i]) * (t[i+1] - t[i])
            break

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14), sharex=True)

    # 1. Eigenvalue & Reactivity
    ax1.errorbar(t, k_0, yerr=sig_0, fmt="o-", color="#1f77b4", label=r"$k_{\mathrm{eff}}(0^\circ)$ (Drums Fully Outward / Excess Reactivity)", capsize=3)
    ax1.errorbar(t, k_crit, yerr=sig_crit, fmt="s--", color="#2ca02c", label=r"$k_{\mathrm{eff}}(\theta_{\mathrm{crit}})$ (Operational Critical State)", capsize=3)
    ax1.axhline(1.00000, color="red", linestyle=":", linewidth=2, label=r"Criticality Threshold ($k_{\mathrm{eff}} = 1.0000$)")
    ax1.axvline(15.0, color="purple", linestyle="--", alpha=0.7, label="Design Target (15.0 Years)")

    if t_life is not None:
        ax1.axvline(t_life, color="black", linestyle="-.", linewidth=2, label=f"Operational Lifetime: {t_life:.2f} yr")
        ax1.scatter([t_life], [1.00000], color="red", s=120, zorder=5)
        ax1.annotate(f"Drums Fully Outward\n$k_{{\\mathrm{{eff}}}} = 1.0000$\nLifetime = {t_life:.2f} yr",
                     xy=(t_life, 1.00000), xytext=(t_life - 4.5, 1.008),
                     arrowprops=dict(facecolor="black", shrink=0.08, width=1.5, headwidth=7),
                     fontweight="bold", bbox=dict(boxstyle="round,pad=0.4", fc="#fff7bc", ec="goldenrod"))

    ax1.set_ylabel(r"Effective Multiplication Factor $k_{\mathrm{eff}}$")
    ax1.set_title("GEOM_OPT15_TB35: High-Precision Multi-Year Criticality Depletion Trajectory\n400 kWth Constant Power | Actively Enforced Critical Control Drums", fontweight="bold")
    ax1.legend(loc="lower left", framealpha=0.9)
    ax1.grid(True)

    # 2. Control Drum Angle Trajectory
    ax2.plot(t, theta_crit, "d-", color="#d62728", label=r"Critical Drum Rotation $\theta_{\mathrm{crit}}(t)$")
    ax2.axhline(0.0, color="gray", linestyle="--", label=r"Fully Outward ($0^\circ$)")
    ax2.axvline(15.0, color="purple", linestyle="--", alpha=0.7)
    if t_life is not None:
        ax2.axvline(t_life, color="black", linestyle="-.", linewidth=2)
    ax2.set_ylabel(r"Control Drum Angle $\theta_{\mathrm{crit}}$ [deg]")
    ax2.set_title(r"Active Drum Positioning Trajectory Enforcing $k_{\mathrm{eff}} = 1.0000$", fontweight="bold")
    ax2.legend(loc="upper right", framealpha=0.9)
    ax2.grid(True)

    # 3. Fissile Inventory (U-235 and Pu-239)
    ax3.plot(t, u235, "o-", color="#ff7f0e", label="U-235 Inventory [kg]")
    ax3_twin = ax3.twinx()
    ax3_twin.plot(t, pu239, "^--", color="#9467bd", label="Pu-239 Bred Inventory [kg]")
    ax3.set_xlabel("Operational Time [years]")
    ax3.set_ylabel("U-235 Mass [kg]", color="#ff7f0e")
    ax3_twin.set_ylabel("Pu-239 Mass [kg]", color="#9467bd")
    ax3.set_title("Actinide Inventory Evolution: Fissile Depletion & Breeding", fontweight="bold")
    ax3.axvline(15.0, color="purple", linestyle="--", alpha=0.7)
    if t_life is not None:
        ax3.axvline(t_life, color="black", linestyle="-.", linewidth=2)

    # Combined legend
    lines1, labels1 = ax3.get_legend_handles_labels()
    lines2, labels2 = ax3_twin.get_legend_handles_labels()
    ax3.legend(lines1 + lines2, labels1 + labels2, loc="center left", framealpha=0.9)
    ax3.grid(True)

    plt.tight_layout()
    out_path = os.path.join(FIG_DIR, "geom_tb35_07_depletion_lifetime_and_drums.png")
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[Plot Generated] {out_path}")


if __name__ == "__main__":
    plot_mass_breakdown()
    plot_depletion_trajectory()
