#!/usr/bin/env python3
"""
generate_neutronic_figures_tb35.py

Generates publication-quality technical figures (300 DPI) for:
  - geom_tb35_05_control_drum_calibration_curve.png: S-curve and differential drum worth
  - geom_tb35_06_reactivity_feedbacks_breakdown.png: Full breakdown of reactivity feedback coefficients
"""

import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/best_configurations/GEOM_OPT15_TB35"
DATA_DIR = os.path.join(BASE_DIR, "data")
FIG_DIR = os.path.join(BASE_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

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

# 1. Drum Calibration Figure
calib_csv = os.path.join(DATA_DIR, "tb35_control_drum_calibration.csv")
df_calib = pd.read_csv(calib_csv)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5), dpi=300)

# Panel 1: k_eff vs Theta
angles = df_calib["theta_deg"].values
keffs = df_calib["k_eff"].values
sigmas = df_calib["sigma_k"].values

ax1.errorbar(angles, keffs, yerr=sigmas, fmt="o-", color="#1b9e77", ecolor="black",
             elinewidth=1.5, capsize=4, markersize=7, linewidth=2.0, label=r"$k_{\mathrm{eff}}(\theta)$ (All 10 Drums)")
ax1.axhline(1.000, color="black", linestyle="--", linewidth=1.5, label=r"Criticality Threshold ($k=1.0000$)")
ax1.plot(80.50, 0.99844, 's', color="#d95f02", markersize=9, label=r"BOL Critical Angle: $\theta_{\mathrm{crit}}=80.50^\circ$ ($k=0.99844$)")
ax1.plot(180.0, 0.97521, '^', color="#7570b3", markersize=9, label=r"Stuck Drum: $k_{\mathrm{stuck}}=0.97521$ ($\mathrm{SDM}=2542$ pcm)")

ax1.set_xlabel(r"Control Drum Rotation Angle $\theta$ [deg]")
ax1.set_ylabel(r"Effective Multiplication Factor $k_{\mathrm{eff}}$")
ax1.set_title(r"(a) Multiplication Factor vs Drum Rotation Angle")
ax1.set_xlim(-5, 185)
ax1.set_ylim(0.965, 1.025)
ax1.legend(loc="upper right", frameon=True, fontsize=9)

for th, k in zip([0.0, 60.0, 120.0, 180.0], [keffs[0], keffs[2], keffs[4], keffs[6]]):
    ax1.annotate(f"{k:.5f}", (th, k), textcoords="offset points", xytext=(0, 10 if th < 90 else -16),
                 ha="center", fontsize=9, fontweight="bold")

# Panel 2: Reactivity Inserted & Differential Worth
worth = df_calib["worth_inserted_pcm"].values

# Spline / smooth curve
th_smooth = np.linspace(0, 180, 200)
# S-curve: W_tot * (th/180 - sin(2*pi*th/180)/(2*pi))
W_tot = worth[-1]
worth_smooth = W_tot * (th_smooth / 180.0 - np.sin(2.0 * math.pi * th_smooth / 180.0) / (2.0 * math.pi))
diff_worth_smooth = (W_tot / 180.0) * (1.0 - np.cos(2.0 * math.pi * th_smooth / 180.0)) # pcm/deg

ax2.plot(angles, worth, 'o', color="#e7298a", markersize=7, label="Monte Carlo Integral Worth [pcm]")
ax2.plot(th_smooth, worth_smooth, '-', color="#e7298a", linewidth=2.0, label="Integral Calibration S-Curve")
ax2.set_xlabel(r"Control Drum Rotation Angle $\theta$ [deg]")
ax2.set_ylabel(r"Integral Absorber Worth Inserted $\Delta\rho$ [pcm]", color="#e7298a")
ax2.tick_params(axis='y', labelcolor="#e7298a")
ax2.set_xlim(-5, 185)
ax2.set_ylim(-200, 5200)

ax2_twin = ax2.twinx()
ax2_twin.plot(th_smooth, diff_worth_smooth, '--', color="#386cb0", linewidth=1.8, label=r"Differential Worth $d\rho/d\theta$ [pcm/deg]")
ax2_twin.set_ylabel(r"Differential Drum Worth [pcm/deg]", color="#386cb0")
ax2_twin.tick_params(axis='y', labelcolor="#386cb0")
ax2_twin.set_ylim(0, 60)

ax2.set_title(r"(b) Integral and Differential Control Drum Worth")
lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax2_twin.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2, loc="center left", frameon=True, fontsize=9)

plt.tight_layout()
fig5_path = os.path.join(FIG_DIR, "geom_tb35_05_control_drum_calibration_curve.png")
plt.savefig(fig5_path, dpi=300)
plt.close()
print(f"[+] Saved: {fig5_path}")

# 2. Reactivity Feedback Breakdown Figure
fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

coeff_labels = [
    r"Fuel Doppler $\alpha_{\mathrm{fuel}}$",
    r"Moderator Total $\alpha_{\mathrm{mod}}$",
    r"Cladding $\alpha_{\mathrm{clad}}$",
    r"Reflector $\alpha_{\mathrm{refl}}$",
    r"TOTAL Isothermal $\alpha_{\mathrm{iso}}$"
]

coeff_vals = [-0.9880, -1.1294, -0.6686, 0.0119, -2.6894]
coeff_errs = [0.1626, 0.3229, 0.2909, 0.2872, 0.5497]

y_pos = np.arange(len(coeff_labels))
colors = ["#d95f02", "#7570b3", "#66a61e", "#e6ab02", "#e41a1c"]

bars = ax.barh(y_pos, coeff_vals, xerr=coeff_errs, color=colors, height=0.55, edgecolor="black", linewidth=1.2,
               error_kw=dict(ecolor="black", lw=1.5, capsize=5, capthick=1.5))

ax.axvline(0.0, color="black", linestyle="--", linewidth=1.2)
ax.set_yticks(y_pos)
ax.set_yticklabels(coeff_labels, fontsize=11, fontweight="bold")
ax.invert_yaxis()  # Top-down
ax.set_xlabel(r"Reactivity Feedback Coefficient [pcm / K]", fontsize=12)
ax.set_title("GEOM_OPT15_TB35: Temperature Reactivity Feedback Breakdown\n"
             r"Total Isothermal Coefficient $\alpha_{\mathrm{iso}} = -2.689 \pm 0.550$ pcm/K (Passive Safety Confirmed)",
             fontsize=12, fontweight="bold")
ax.set_xlim(-3.6, 0.8)

for bar, val, err in zip(bars, coeff_vals, coeff_errs):
    x_pos = val - err - 0.12 if val < 0 else val + err + 0.05
    ha = "right" if val < 0 else "left"
    ax.text(x_pos, bar.get_y() + bar.get_height()/2.0, f"{val:+.3f} ± {err:.3f} pcm/K",
            va="center", ha=ha, fontsize=10, fontweight="bold")

# Inset annotation
ax.text(-3.4, 4.3, r"Prompt Power Coefficient: $\alpha_P = -0.744 \pm 0.116$ pcm/kWth" + "\n" +
                   r"Moderator Density Coeff: $\alpha_{\mathrm{mod,dens}} = +310.9 \pm 14.4$ pcm/%",
        fontsize=10, bbox=dict(boxstyle="round,pad=0.5", facecolor="#fff9c4", edgecolor="#fbc02d", lw=1.2))

plt.tight_layout()
fig6_path = os.path.join(FIG_DIR, "geom_tb35_06_reactivity_feedbacks_breakdown.png")
plt.savefig(fig6_path, dpi=300)
plt.close()
print(f"[+] Saved: {fig6_path}")
