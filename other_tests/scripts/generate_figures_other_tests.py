#!/usr/bin/env python3
"""
generate_figures_other_tests.py

Generates publication-quality technical figures for the Moderator Thermal Buffer study
in 'other_tests'.
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/other_tests"
DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

csv_path = os.path.join(DATA_DIR, "parametric_campaign_results.csv")
df = pd.read_csv(csv_path)

plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.titlesize": 14
})

# -----------------------------------------------------------------------------
# FIGURE 1: Thermal Buffer & Moderator Operating Temperature
# -----------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)

# Panel 1: Peak Moderator Temperature across Cases
cases_labels = [
    "Case 1\nBaseline\n(920K, 1.5mm)",
    "Case 2\nWide Holes\n(920K, 1.8mm)",
    "Case 3\nLower Inlet\n(900K, 1.5mm)",
    "Case 3B\nTarget Buffer\n(892K, 1.5mm)",
    "Case 4\nCombined\n(905K, 1.7mm)",
    "Case 5\nHigh-Buffer\n(895K, 1.8mm)"
]

# Re-order to match logical progression: 1, 2, 3, 3B, 4, 5
order = ["case1_baseline", "case2_wide_channel", "case3_inlet_lowered", "case3b_opt_target_892K", "case4_synergistic_balanced", "case5_target_high_buffer"]
sub_df = df.set_index("case_id").loc[order].reset_index()

t_mod = sub_df["T_mod_peak_K"].values
buffers = sub_df["buffer_mod_K"].values

colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

bars1 = ax1.bar(range(len(cases_labels)), t_mod, color=colors, width=0.55, edgecolor="black", linewidth=1.2)
ax1.axhline(950.0, color="crimson", linestyle="--", linewidth=2.0, label="Dissociation Limit (950.0 K)")
ax1.axhline(920.0, color="gray", linestyle=":", linewidth=1.5, label="Target Ceiling (920.0 K)")
ax1.axhline(915.0, color="green", linestyle=":", linewidth=1.5, label="Target Floor (915.0 K)")
ax1.set_xticks(range(len(cases_labels)))
ax1.set_xticklabels(cases_labels, rotation=0)
ax1.set_ylabel("Peak Canned YH2 Temperature [K]")
ax1.set_title("Peak Moderator Temperature Comparison")
ax1.set_ylim(890.0, 960.0)
ax1.legend(loc="upper right", frameon=True)

for bar in bars1:
    h = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2.0, h + 1.0, f"{h:.1f} K", ha="center", va="bottom", fontweight="bold", fontsize=9.5)

# Panel 2: Safety Buffer to 950 K Limit
bars2 = ax2.bar(range(len(cases_labels)), buffers, color=colors, width=0.55, edgecolor="black", linewidth=1.2)
ax2.axhline(35.0, color="crimson", linestyle="--", linewidth=2.0, label="Target Requirement (+35.0 K)")
ax2.axhline(10.0, color="orange", linestyle=":", linewidth=1.5, label="Baseline Margin (+10.0 K)")
ax2.set_xticks(range(len(cases_labels)))
ax2.set_xticklabels(cases_labels, rotation=0)
ax2.set_ylabel("Thermal Buffer to 950 K Ceiling [K]")
ax2.set_title("Thermal Safety Margin (+Delta T to Dissociation)")
ax2.set_ylim(0.0, 45.0)
ax2.legend(loc="upper left", frameon=True)

for bar in bars2:
    h = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2.0, h + 0.8, f"+{h:.1f} K", ha="center", va="bottom", fontweight="bold", fontsize=9.5)

plt.tight_layout()
fig1_path = os.path.join(FIG_DIR, "fig1_thermal_buffer_vs_inlet_and_channel_radius.png")
plt.savefig(fig1_path, dpi=300)
plt.close()
print(f"[+] Saved: {fig1_path}")

# -----------------------------------------------------------------------------
# FIGURE 2: Neutronic Impact of Channel Widening vs Inlet Lowering
# -----------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)

# Panel 1: k_eff(0 deg) vs r_add
r_add_vals = [1.50, 1.70, 1.80]
keff_r = [
    sub_df.loc[sub_df["case_id"] == "case1_baseline", "k_eff_0deg"].values[0],
    sub_df.loc[sub_df["case_id"] == "case4_synergistic_balanced", "k_eff_0deg"].values[0],
    sub_df.loc[sub_df["case_id"] == "case2_wide_channel", "k_eff_0deg"].values[0]
]
err_r = [
    sub_df.loc[sub_df["case_id"] == "case1_baseline", "sigma_k_0deg"].values[0],
    sub_df.loc[sub_df["case_id"] == "case4_synergistic_balanced", "sigma_k_0deg"].values[0],
    sub_df.loc[sub_df["case_id"] == "case2_wide_channel", "sigma_k_0deg"].values[0]
]

ax1.errorbar(r_add_vals, keff_r, yerr=err_r, fmt="o-", color="#d62728", ecolor="black", elinewidth=1.5, capsize=5, markersize=8, linewidth=2.0, label="k_eff(0 deg) at BOL")
ax1.axhline(1.000, color="black", linestyle="--", linewidth=1.5, label="Criticality Threshold (k = 1.000)")
ax1.axhline(1.016, color="blue", linestyle=":", linewidth=1.5, label="15-Year Depletion Requirement (k >= 1.016)")
ax1.set_xlabel("Auxiliary Channel Radius r_add [mm]")
ax1.set_ylabel("Multiplication Factor k_eff (0 deg)")
ax1.set_title("Reactivity Loss from Moderator Voiding (YH2 Mass Trimming)")
ax1.set_ylim(0.985, 1.025)
ax1.legend(loc="lower left", frameon=True)

for x, y in zip(r_add_vals, keff_r):
    ax1.annotate(f"k = {y:.5f}", (x, y), textcoords="offset points", xytext=(0, 10), ha="center", fontweight="bold", fontsize=9.5)

# Panel 2: Projected Operational Lifespan across Cases
lifetimes = sub_df["lifetime_years"].values
bar_colors = ["#2ca02c" if l >= 15.0 else "#d62728" for l in lifetimes]

bars_life = ax2.bar(range(len(cases_labels)), [max(0.0, l) for l in lifetimes], color=bar_colors, width=0.55, edgecolor="black", linewidth=1.2)
ax2.axhline(15.0, color="black", linestyle="--", linewidth=2.0, label="15.0-Year Mission Requirement")
ax2.set_xticks(range(len(cases_labels)))
ax2.set_xticklabels(cases_labels, rotation=0)
ax2.set_ylabel("Projected Core Lifespan [Years]")
ax2.set_title("Operational Lifetime Comparison (Depletion Reserve)")
ax2.set_ylim(0.0, 22.0)
ax2.legend(loc="upper right", frameon=True)

for idx, bar in enumerate(bars_life):
    l = lifetimes[idx]
    if l >= 0:
        ax2.text(bar.get_x() + bar.get_width()/2.0, bar.get_height() + 0.5, f"{l:.1f} yr", ha="center", va="bottom", fontweight="bold", fontsize=9.5)
    else:
        ax2.text(bar.get_x() + bar.get_width()/2.0, 0.5, "SUBCRITICAL\n(0 yr)", ha="center", va="bottom", fontweight="bold", fontsize=8.5, color="crimson")

plt.tight_layout()
fig2_path = os.path.join(FIG_DIR, "fig2_neutronic_penalty_of_channel_widening.png")
plt.savefig(fig2_path, dpi=300)
plt.close()
print(f"[+] Saved: {fig2_path}")

# -----------------------------------------------------------------------------
# FIGURE 3: Flight Configuration Performance Spider / Bar Metric Audit
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), dpi=300)

metrics = [
    ("buffer_mod_K", "Moderator Thermal Buffer [K]", 35.0, "Target: >= +35 K", "#2ca02c", (0, 45)),
    ("p_h2_kpa", "Sieverts H2 Pressure [kPa]", 31.8, "Ref Ceiling: 31.8 kPa", "#1f77b4", (0, 3.0)),
    ("M_tot_kg", "Core Subassembly Mass [kg]", 600.0, "Ceiling: 600.0 kg", "#9467bd", (520, 620)),
    ("dp_core_kPa", "Core Pressure Drop [kPa]", 35.0, "Ceiling: 35.0 kPa", "#ff7f0e", (0, 25)),
    ("eta_brayton_pct", "Brayton Cycle Efficiency [%]", 30.0, "Target: >= 30.0%", "#17becf", (25, 45)),
    ("lifetime_years", "Operational Lifespan [yr]", 15.0, "Target: >= 15.0 yr", "#d62728", (0, 22))
]

comp_cases = ["case1_baseline", "case3_inlet_lowered", "case3b_opt_target_892K"]
comp_labels = ["Baseline (920K)", "Tin=900K", "Case 3B (Tin=892K, +35K)"]
comp_df = df.set_index("case_id").loc[comp_cases].reset_index()

bar_colors_comp = ["#7f7f7f", "#1f77b4", "#2ca02c"]

for ax, (metric, label, limit, lim_label, c, ylim) in zip(axes.flatten(), metrics):
    vals = comp_df[metric].values
    bars = ax.bar(range(len(comp_labels)), vals, color=bar_colors_comp, width=0.55, edgecolor="black", linewidth=1.2)
    ax.axhline(limit, color="crimson", linestyle="--", linewidth=1.8, label=lim_label)
    ax.set_xticks(range(len(comp_labels)))
    ax.set_xticklabels(comp_labels, fontsize=9.5)
    ax.set_ylabel(label)
    ax.set_title(label.split("[")[0].strip())
    ax.set_ylim(ylim)
    ax.legend(loc="best", fontsize=9, frameon=True)
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2.0, h + (ylim[1]-ylim[0])*0.02, f"{h:.2f}", ha="center", va="bottom", fontweight="bold", fontsize=9)

plt.tight_layout()
fig3_path = os.path.join(FIG_DIR, "fig3_flight_configuration_comparison.png")
plt.savefig(fig3_path, dpi=300)
plt.close()
print(f"[+] Saved: {fig3_path}")
