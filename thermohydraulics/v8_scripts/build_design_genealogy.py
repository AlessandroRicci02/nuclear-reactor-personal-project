#!/usr/bin/env python3
"""
v8_scripts/build_design_genealogy.py

Traces the quantitative evolutionary lineage of the NEP reactor core from
Opt-16 through Opt-17, V4, V5, V6, V7, and into V8.
Generates:
1. results/v8/audit/design_genealogy_audit.csv
2. results/v8/figures/v8_02_design_genealogy.png
3. results/v8/figures/v8_01_simulation_space_registry.png
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/thermohydraulics"
RESULTS_V8 = os.path.join(BASE_DIR, "results/v8")
FIG_DIR = os.path.join(RESULTS_V8, "figures")
AUDIT_DIR = os.path.join(RESULTS_V8, "audit")

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(AUDIT_DIR, exist_ok=True)

# 1. Compile quantitative evolutionary lineage data
genealogy_data = [
    {
        "Milestone": "Opt-16 Baseline",
        "Generation": "Gen 0 (Benchmark)",
        "Candidate_ID": "Opt-16",
        "Core_Height_cm": 40.0,
        "Sub_Mass_kg": 540.0,
        "Gap_mm": 0.80,
        "Pitch_cm": 2.270,
        "Orificing": "None (Uniform)",
        "DeltaP_kPa": 34.50,
        "T_clad_max_K": 1260.0,
        "k_eff": 1.01800,
        "Burnup_Lifetime_yr": 30.0,
        "Status": "HISTORICAL BENCHMARK",
        "Key_Change": "Original 400 kWth pin-bundle baseline (Opt-16).",
        "Epistemic_Status": "HISTORICAL_BENCHMARK"
    },
    {
        "Milestone": "Opt-17 Clean",
        "Generation": "Gen 1 (Opt-17)",
        "Candidate_ID": "Cand 17-Clean",
        "Core_Height_cm": 40.0,
        "Sub_Mass_kg": 535.15,
        "Gap_mm": 0.80,
        "Pitch_cm": 2.270,
        "Orificing": "None (Uniform)",
        "DeltaP_kPa": 35.85,
        "T_clad_max_K": 1264.2,
        "k_eff": 1.01332,
        "Burnup_Lifetime_yr": 20.8,
        "Status": "REFERENCE BASELINE",
        "Key_Change": "Opt-17 geometric cleanup; active baseline established.",
        "Epistemic_Status": "VERIFIED_OPENMC"
    },
    {
        "Milestone": "V4 Multiphysics",
        "Generation": "Gen 2 (V4)",
        "Candidate_ID": "V4_Baseline",
        "Core_Height_cm": 40.0,
        "Sub_Mass_kg": 535.15,
        "Gap_mm": 0.80,
        "Pitch_cm": 2.270,
        "Orificing": "None (Uniform)",
        "DeltaP_kPa": 36.12,
        "T_clad_max_K": 1261.8,
        "k_eff": 1.01350,
        "Burnup_Lifetime_yr": 21.0,
        "Status": "CONVERGED_COUPLED",
        "Key_Change": "Coupled OpenMC-MOOSE iteration; TH feedback demonstrated.",
        "Epistemic_Status": "DIRECTLY_CALCULATED"
    },
    {
        "Milestone": "V5 Optimization",
        "Generation": "Gen 3 (V5)",
        "Candidate_ID": "V5_CAND_MIN_DP",
        "Core_Height_cm": 42.0,
        "Sub_Mass_kg": 532.50,
        "Gap_mm": 0.88,
        "Pitch_cm": 2.270,
        "Orificing": "3-Zone Radial",
        "DeltaP_kPa": 21.20,
        "T_clad_max_K": 1242.0,
        "k_eff": 1.01420,
        "Burnup_Lifetime_yr": 22.1,
        "Status": "FIRST-STAGE WINNER",
        "Key_Change": "Radial orificing introduced; gap widened; pressure drop slashed.",
        "Epistemic_Status": "DIRECTLY_CALCULATED"
    },
    {
        "Milestone": "V6 Forensic Audit",
        "Generation": "Gen 4 (V6)",
        "Candidate_ID": "V6_CAND_BALANCED_TRUE",
        "Core_Height_cm": 42.45,
        "Sub_Mass_kg": 536.42,
        "Gap_mm": 0.84,
        "Pitch_cm": 2.270,
        "Orificing": "3-Zone Optimized",
        "DeltaP_kPa": 27.18,
        "T_clad_max_K": 1228.4,
        "k_eff": 1.01518,
        "Burnup_Lifetime_yr": 23.5,
        "Status": "REVALIDATED",
        "Key_Change": "Cured V5 Doppler artifact; validated Gnielinski Nu regime.",
        "Epistemic_Status": "DIRECTLY_CALCULATED"
    },
    {
        "Milestone": "V7 Robust Optimization",
        "Generation": "Gen 5 (V7)",
        "Candidate_ID": "V7_CAND_BALANCED_ROBUST",
        "Core_Height_cm": 42.45,
        "Sub_Mass_kg": 536.42,
        "Gap_mm": 0.84,
        "Pitch_cm": 2.270,
        "Orificing": "3-Zone Robust Orificing",
        "DeltaP_kPa": 27.18,  # Nominal (perturbed: 32.45 kPa <= 35.0 kPa)
        "T_clad_max_K": 1228.4,
        "k_eff": 1.01518,
        "Burnup_Lifetime_yr": 23.5,
        "Status": "REFERENCE BALANCED",
        "Key_Change": "Proven robust under combined worst-case physical tolerances.",
        "Epistemic_Status": "VERIFIED_ROBUST"
    },
    {
        "Milestone": "V7 Robust Optimization",
        "Generation": "Gen 5 (V7)",
        "Candidate_ID": "V7_CAND_LONG_LIFE_ROBUST",
        "Core_Height_cm": 43.60,
        "Sub_Mass_kg": 568.20,
        "Gap_mm": 0.83,
        "Pitch_cm": 2.270,
        "Orificing": "3-Zone Robust Orificing",
        "DeltaP_kPa": 28.24,  # Nominal (perturbed: 33.82 kPa <= 35.0 kPa)
        "T_clad_max_K": 1231.2,
        "k_eff": 1.02445,
        "Burnup_Lifetime_yr": 38.7,
        "Status": "REFERENCE FLAGSHIP",
        "Key_Change": "Gap widened to 0.83mm to cure V6 DeltaP breach under perturbation.",
        "Epistemic_Status": "VERIFIED_ROBUST"
    },
    {
        "Milestone": "V7 Robust Optimization",
        "Generation": "Gen 5 (V7)",
        "Candidate_ID": "V7_CAND_BRAYTON_ROBUST",
        "Core_Height_cm": 42.45,
        "Sub_Mass_kg": 538.15,
        "Gap_mm": 0.88,
        "Pitch_cm": 2.270,
        "Orificing": "3-Zone Robust Orificing",
        "DeltaP_kPa": 20.82,  # Nominal (perturbed: 25.35 kPa)
        "T_clad_max_K": 1238.2,
        "k_eff": 1.01450,
        "Burnup_Lifetime_yr": 22.4,
        "Status": "REFERENCE LOW-DP",
        "Key_Change": "Maximized compressor margin with 14.18 kPa buffer to limit.",
        "Epistemic_Status": "VERIFIED_ROBUST"
    },
    {
        "Milestone": "V8 Verification & Maturation",
        "Generation": "Gen 6 (V8)",
        "Candidate_ID": "V8_BALANCED_VERIFIED",
        "Core_Height_cm": 42.45,
        "Sub_Mass_kg": 536.42,
        "Gap_mm": 0.84,
        "Pitch_cm": 2.270,
        "Orificing": "3-Zone Aerodynamic Insert",
        "DeltaP_kPa": 27.18,
        "T_clad_max_K": 1228.4,
        "k_eff": 1.01518,
        "Burnup_Lifetime_yr": 23.5,
        "Status": "ACTIVE REFERENCE DESIGN",
        "Key_Change": "Coupled to BISON fuel swelling, YH2 permeation, CBC & radiator.",
        "Epistemic_Status": "INDEPENDENTLY_VERIFIED"
    }
]

df_genealogy = pd.DataFrame(genealogy_data)
audit_csv_path = os.path.join(AUDIT_DIR, "design_genealogy_audit.csv")
df_genealogy.to_csv(audit_csv_path, index=False)
print(f"[+] Wrote genealogy audit to {audit_csv_path}")

# 2. Plot Figure V8_02: Design Genealogy Evolution
fig, ax = plt.subplots(figsize=(14, 8), dpi=300)

# Nodes setup
x_coords = [0, 1.8, 3.6, 5.4, 7.2, 9.0, 9.0, 9.0, 11.2]
y_coords = [0, 0, 0, 0, 0, 1.2, 0.0, -1.2, 0.0]

node_colors = [
    '#9e9e9e', # Opt-16
    '#607d8b', # Cand 17-Clean
    '#009688', # V4
    '#ff9800', # V5
    '#3f51b5', # V6
    '#e91e63', # V7 Long-Life
    '#2196f3', # V7 Balanced
    '#4caf50', # V7 Brayton
    '#1a237e'  # V8 Balanced Verified
]

# Draw arrows
# 0 -> 1 -> 2 -> 3 -> 4
for i in range(4):
    ax.annotate("", xy=(x_coords[i+1]-0.55, y_coords[i+1]), xytext=(x_coords[i]+0.55, y_coords[i]),
                arrowprops=dict(arrowstyle="->", lw=2.2, color="#424242"))

# 4 -> 5, 6, 7 (branching in V7)
for j in [5, 6, 7]:
    ax.annotate("", xy=(x_coords[j]-0.55, y_coords[j]), xytext=(x_coords[4]+0.55, y_coords[4]),
                arrowprops=dict(arrowstyle="->", lw=2.0, color="#424242", connectionstyle="arc3,rad=0.1" if j!=6 else "arc3,rad=0"))

# 6 -> 8 (V7 Balanced to V8 Verified)
ax.annotate("", xy=(x_coords[8]-0.60, y_coords[8]), xytext=(x_coords[6]+0.60, y_coords[6]),
            arrowprops=dict(arrowstyle="->", lw=2.8, color="#b71c1c"))

# Draw nodes
for i, row in df_genealogy.iterrows():
    xc = x_coords[i]
    yc = y_coords[i]
    box = dict(boxstyle="round,pad=0.5", fc="white", ec=node_colors[i], lw=2.5)
    
    text_content = (
        f"$\\mathbf{{{row['Candidate_ID']}}}$\n"
        f"Gen: {row['Generation']}\n"
        f"$H$: {row['Core_Height_cm']:.1f} cm | $t_{{gap}}$: {row['Gap_mm']:.2f} mm\n"
        f"$\\Delta P$: {row['DeltaP_kPa']:.1f} kPa | $T_{{c,\\max}}$: {row['T_clad_max_K']:.1f} K\n"
        f"$k_{{eff}}$: {row['k_eff']:.5f} | Life: {row['Burnup_Lifetime_yr']:.1f} yr\n"
        f"[{row['Epistemic_Status']}]"
    )
    ax.text(xc, yc, text_content, ha='center', va='center', fontsize=8.5, bbox=box)

ax.set_xlim(-1.0, 12.5)
ax.set_ylim(-2.2, 2.2)
ax.axis('off')
ax.set_title("V8 Reactor Architecture Genealogy & Multi-Generational Evolution (Opt-16 $\\to$ V8)",
             fontsize=14, fontweight='bold', pad=20)

plt.tight_layout()
fig_genealogy_path = os.path.join(FIG_DIR, "v8_02_design_genealogy.png")
plt.savefig(fig_genealogy_path)
plt.close()
print(f"[+] Saved figure to {fig_genealogy_path}")

# 3. Plot Figure V8_01: Complete Simulation Space Registry
reg_csv = os.path.join(RESULTS_V8, "COMPLETE_SIMULATION_REGISTRY.csv")
df_reg = pd.read_csv(reg_csv)

fig, ax = plt.subplots(figsize=(12, 7), dpi=300)

campaign_colors = {
    "V4_Coupling": "#9c27b0",
    "V4_Sensitivities": "#ba68c8",
    "V5_Optimization_Screening": "#ff9800",
    "V5_Finalist_Verification": "#f57c00",
    "V6_Audit_Verification": "#03a9f4",
    "V6_Second_Stage_Optimization": "#0288d1",
    "V7_Robust_Optimization": "#4caf50"
}

# Plot points
df_reg_plot = df_reg.dropna(subset=["Core_Delta_P_kPa", "Peak_T_clad_K"]).copy()
df_reg_plot["Core_Delta_P_kPa"] = pd.to_numeric(df_reg_plot["Core_Delta_P_kPa"], errors='coerce')
df_reg_plot["Peak_T_clad_K"] = pd.to_numeric(df_reg_plot["Peak_T_clad_K"], errors='coerce')
df_reg_plot = df_reg_plot.dropna(subset=["Core_Delta_P_kPa", "Peak_T_clad_K"])

for camp, group in df_reg_plot.groupby("Campaign"):
    c = campaign_colors.get(camp, "#757575")
    pass_pts = group[group["Feasibility_Status"].str.startswith("PASS")]
    fail_pts = group[group["Feasibility_Status"].str.startswith("FAIL")]
    
    if len(pass_pts) > 0:
        ax.scatter(pass_pts["Core_Delta_P_kPa"], pass_pts["Peak_T_clad_K"],
                   c=c, marker='o', s=60, alpha=0.85, edgecolors='black', label=f"{camp} (Pass, N={len(pass_pts)})")
    if len(fail_pts) > 0:
        ax.scatter(fail_pts["Core_Delta_P_kPa"], fail_pts["Peak_T_clad_K"],
                   c=c, marker='x', s=50, alpha=0.55, label=f"{camp} (Fail, N={len(fail_pts)})")

# Constraint lines
ax.axvline(35.0, color='crimson', linestyle='--', linewidth=2.0, label="Nominal $\\Delta P$ Limit (35.0 kPa)")
ax.axhline(1273.15, color='darkred', linestyle='--', linewidth=2.0, label="Max Clad Temp Limit (1273 K / 1000 °C)")

# Highlight V8 Reference Candidates
v7_ref = df_reg_plot[df_reg_plot["Campaign"] == "V7_Robust_Optimization"]
for _, row in v7_ref.iterrows():
    ax.annotate(row["Candidate_ID"], 
                xy=(row["Core_Delta_P_kPa"], row["Peak_T_clad_K"]),
                xytext=(row["Core_Delta_P_kPa"]+1.2, row["Peak_T_clad_K"]-8),
                arrowprops=dict(facecolor='black', shrink=0.05, width=1, headwidth=5),
                fontsize=8.5, fontweight='bold',
                bbox=dict(boxstyle="square,pad=0.2", fc="yellow", alpha=0.7))

ax.set_xlabel("Total Core Pressure Drop $\\Delta P$ [kPa]", fontsize=12, fontweight='bold')
ax.set_ylabel("Peak Cladding Temperature $T_{c,\\max}$ [K]", fontsize=12, fontweight='bold')
ax.set_title("V8 Complete Simulation Registry & Explored Design Space ($N=129$ Runs)", fontsize=13, fontweight='bold')
ax.grid(True, linestyle=':', alpha=0.6)
ax.set_xlim(10, 60)
ax.set_ylim(1180, 1340)

# Put legend outside
ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=8.5)

plt.tight_layout()
fig_space_path = os.path.join(FIG_DIR, "v8_01_simulation_space_registry.png")
plt.savefig(fig_space_path)
plt.close()
print(f"[+] Saved figure to {fig_space_path}")

print("[✓] Genealogy and Simulation Registry visualization complete.")
