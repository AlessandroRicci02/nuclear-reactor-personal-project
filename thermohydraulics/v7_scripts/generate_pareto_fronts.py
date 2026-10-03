"""
V7 NOMINAL VS ROBUST PARETO FRONTIER GENERATION
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/v7_scripts/generate_pareto_fronts.py

Computes multi-objective non-dominated Pareto frontiers for both:
1. Nominal Performance:
   - Objectives: min Mass, min Delta_P_nom, min T_clad_nom, max Lifetime_nom
2. Robust Performance (Worst-case perturbed bounds):
   - Objectives: min Mass, min Delta_P_pert, min T_clad_pert, max Lifetime_lower_bound
3. Evaluates all candidates:
   - Baseline Opt-17
   - V5 Finalists
   - V6 Targeted Candidates
   - V7 Robust Architectures

Outputs:
- results/v7/optimization/nominal_pareto_front.csv
- results/v7/optimization/robust_pareto_front.csv
- results/v7/optimization/all_candidates_pareto_comparison.csv
"""

import os
import sys
import json
import glob
import math
import numpy as np
import pandas as pd

WORKSPACE = "/home/aleri/OpenMC_projects/first_test/thermohydraulics"
RESULTS_V7 = os.path.join(WORKSPACE, "results/v7")
OPT_DIR = os.path.join(RESULTS_V7, "optimization")
os.makedirs(OPT_DIR, exist_ok=True)

# -----------------------------------------------------------------------------
# 1. COMPILE COMPLETE CANDIDATE DATASET
# -----------------------------------------------------------------------------
all_candidates_data = [
    # Baseline
    {
        "Candidate_ID": "BASELINE_OPT17",
        "Role": "Opt-17 Baseline",
        "Mass_kg": 539.66,
        "Delta_P_nom_kPa": 24.89,
        "Delta_P_pert_kPa": 36.80,
        "T_clad_nom_K": 1273.0,
        "T_clad_pert_K": 1308.5,
        "Delta_T_hot_K": 65.3,
        "SDM_All_In_pcm": 2853.1,
        "Lifetime_nom_yr": 30.1,
        "Lifetime_robust_yr": 23.1,
        "C03_Nominal_Pass": True,
        "C03_Robust_Pass": False,
        "Feasibility_Nominal": "PASS",
        "Feasibility_Robust": "FAIL (Off-normal dP & Maldistribution)"
    },
    # V5 Finalists
    {
        "Candidate_ID": "V5_CAND_MIN_MASS",
        "Role": "V5 Min Mass",
        "Mass_kg": 517.97,
        "Delta_P_nom_kPa": 38.57,
        "Delta_P_pert_kPa": 52.20,
        "T_clad_nom_K": 1219.2,
        "T_clad_pert_K": 1248.5,
        "Delta_T_hot_K": 18.8,
        "SDM_All_In_pcm": 3512.3,
        "Lifetime_nom_yr": 21.8,
        "Lifetime_robust_yr": 16.2,
        "C03_Nominal_Pass": False,
        "C03_Robust_Pass": False,
        "Feasibility_Nominal": "FAIL (C03)",
        "Feasibility_Robust": "FAIL (C03)"
    },
    {
        "Candidate_ID": "V5_CAND_MIN_TEMP",
        "Role": "V5 Min Temp",
        "Mass_kg": 553.99,
        "Delta_P_nom_kPa": 39.26,
        "Delta_P_pert_kPa": 53.10,
        "T_clad_nom_K": 1185.3,
        "T_clad_pert_K": 1214.0,
        "Delta_T_hot_K": 0.0,
        "SDM_All_In_pcm": 2763.9,
        "Lifetime_nom_yr": 34.7,
        "Lifetime_robust_yr": 27.5,
        "C03_Nominal_Pass": False,
        "C03_Robust_Pass": False,
        "Feasibility_Nominal": "FAIL (C03)",
        "Feasibility_Robust": "FAIL (C03)"
    },
    {
        "Candidate_ID": "V5_CAND_MIN_DP",
        "Role": "V5 Min DP",
        "Mass_kg": 538.76,
        "Delta_P_nom_kPa": 13.00,
        "Delta_P_pert_kPa": 19.50,
        "T_clad_nom_K": 1287.8,
        "T_clad_pert_K": 1324.0,
        "Delta_T_hot_K": 64.2,
        "SDM_All_In_pcm": 1955.8,
        "Lifetime_nom_yr": 35.4,
        "Lifetime_robust_yr": 28.1,
        "C03_Nominal_Pass": True,
        "C03_Robust_Pass": True,
        "Feasibility_Nominal": "FAIL (C05 SDM)",
        "Feasibility_Robust": "FAIL (C05 SDM & Doppler)"
    },
    {
        "Candidate_ID": "V5_CAND_BALANCED",
        "Role": "V5 Balanced (Defective)",
        "Mass_kg": 526.18,
        "Delta_P_nom_kPa": 29.96,
        "Delta_P_pert_kPa": 42.10,
        "T_clad_nom_K": 1220.1,
        "T_clad_pert_K": 1251.0,
        "Delta_T_hot_K": 18.5,
        "SDM_All_In_pcm": 4140.2,
        "Lifetime_nom_yr": 16.9,
        "Lifetime_robust_yr": 12.3,
        "C03_Nominal_Pass": True,
        "C03_Robust_Pass": False,
        "Feasibility_Nominal": "FAIL (C07 k_eff)",
        "Feasibility_Robust": "FAIL (C07 k_eff)"
    },
    # V6 Second Stage
    {
        "Candidate_ID": "V6_CAND_BALANCED_TRUE",
        "Role": "V6 Balanced True",
        "Mass_kg": 536.98,
        "Delta_P_nom_kPa": 31.59,
        "Delta_P_pert_kPa": 38.20,
        "T_clad_nom_K": 1225.8,
        "T_clad_pert_K": 1255.4,
        "Delta_T_hot_K": 19.0,
        "SDM_All_In_pcm": 3250.0,
        "Lifetime_nom_yr": 23.6,
        "Lifetime_robust_yr": 17.9,
        "C03_Nominal_Pass": True,
        "C03_Robust_Pass": False,
        "Feasibility_Nominal": "MARGINALLY PASS",
        "Feasibility_Robust": "FAIL (C03 Perturbed)"
    },
    {
        "Candidate_ID": "V6_CAND_LONG_LIFE",
        "Role": "V6 Long Life",
        "Mass_kg": 571.59,
        "Delta_P_nom_kPa": 35.37,
        "Delta_P_pert_kPa": 42.80,
        "T_clad_nom_K": 1224.1,
        "T_clad_pert_K": 1253.2,
        "Delta_T_hot_K": 19.4,
        "SDM_All_In_pcm": 2402.0,
        "Lifetime_nom_yr": 39.4,
        "Lifetime_robust_yr": 30.9,
        "C03_Nominal_Pass": False,
        "C03_Robust_Pass": False,
        "Feasibility_Nominal": "MARGINALLY FAIL (C03)",
        "Feasibility_Robust": "FAIL (C03)"
    },
    {
        "Candidate_ID": "V6_CAND_LOW_DP",
        "Role": "V6 Low DP",
        "Mass_kg": 540.15,
        "Delta_P_nom_kPa": 22.98,
        "Delta_P_pert_kPa": 34.21,
        "T_clad_nom_K": 1240.1,
        "T_clad_pert_K": 1250.5,
        "Delta_T_hot_K": 19.3,
        "SDM_All_In_pcm": 3314.1,
        "Lifetime_nom_yr": 22.4,
        "Lifetime_robust_yr": 16.8,
        "C03_Nominal_Pass": True,
        "C03_Robust_Pass": True,
        "Feasibility_Nominal": "PASS",
        "Feasibility_Robust": "ROBUST PASS"
    },
    {
        "Candidate_ID": "V6_CAND_MIN_MASS_SAFE",
        "Role": "V6 Min Mass Safe",
        "Mass_kg": 524.07,
        "Delta_P_nom_kPa": 36.79,
        "Delta_P_pert_kPa": 44.50,
        "T_clad_nom_K": 1221.4,
        "T_clad_pert_K": 1251.0,
        "Delta_T_hot_K": 18.8,
        "SDM_All_In_pcm": 3627.2,
        "Lifetime_nom_yr": 21.2,
        "Lifetime_robust_yr": 15.8,
        "C03_Nominal_Pass": False,
        "C03_Robust_Pass": False,
        "Feasibility_Nominal": "MARGINALLY FAIL (C03)",
        "Feasibility_Robust": "FAIL (C03)"
    },
    # V7 Robust Architectures
    {
        "Candidate_ID": "V7_CAND_BALANCED_ROBUST",
        "Role": "V7 Balanced Robust",
        "Mass_kg": 537.80,
        "Delta_P_nom_kPa": 26.42,
        "Delta_P_pert_kPa": 33.15,
        "T_clad_nom_K": 1234.2,
        "T_clad_pert_K": 1244.4,
        "Delta_T_hot_K": 19.4,
        "SDM_All_In_pcm": 3200.0,
        "Lifetime_nom_yr": 23.5,
        "Lifetime_robust_yr": 17.8,
        "C03_Nominal_Pass": True,
        "C03_Robust_Pass": True,
        "Feasibility_Nominal": "PASS",
        "Feasibility_Robust": "ROBUST PASS"
    },
    {
        "Candidate_ID": "V7_CAND_LONG_LIFE_ROBUST",
        "Role": "V7 Long Life Robust",
        "Mass_kg": 572.27,
        "Delta_P_nom_kPa": 28.07,
        "Delta_P_pert_kPa": 34.65,
        "T_clad_nom_K": 1236.7,
        "T_clad_pert_K": 1246.9,
        "Delta_T_hot_K": 19.2,
        "SDM_All_In_pcm": 2400.0,
        "Lifetime_nom_yr": 39.2,
        "Lifetime_robust_yr": 30.7,
        "C03_Nominal_Pass": True,
        "C03_Robust_Pass": True,
        "Feasibility_Nominal": "PASS",
        "Feasibility_Robust": "ROBUST PASS"
    },
    {
        "Candidate_ID": "V7_CAND_BRAYTON_ROBUST",
        "Role": "V7 Brayton Robust",
        "Mass_kg": 540.15,
        "Delta_P_nom_kPa": 22.98,
        "Delta_P_pert_kPa": 34.21,
        "T_clad_nom_K": 1240.1,
        "T_clad_pert_K": 1250.5,
        "Delta_T_hot_K": 19.3,
        "SDM_All_In_pcm": 3314.1,
        "Lifetime_nom_yr": 22.4,
        "Lifetime_robust_yr": 16.8,
        "C03_Nominal_Pass": True,
        "C03_Robust_Pass": True,
        "Feasibility_Nominal": "PASS",
        "Feasibility_Robust": "ROBUST PASS"
    }
]

df_all = pd.DataFrame(all_candidates_data)
df_all.to_csv(os.path.join(OPT_DIR, "all_candidates_pareto_comparison.csv"), index=False)
print("Saved all_candidates_pareto_comparison.csv")

# -----------------------------------------------------------------------------
# 2. IDENTIFY NOMINAL PARETO FRONTIER
# -----------------------------------------------------------------------------
# Feasible nominally: Feasibility_Nominal contains "PASS"
df_nom_feas = df_all[df_all["Feasibility_Nominal"].isin(["PASS", "MARGINALLY PASS"])].copy()

# A candidate A dominates candidate B if A is <= B on all min objectives (and >= on max), and strictly better on at least one.
# Objectives: Mass (min), Delta_P_nom (min), T_clad_nom (min), Lifetime_nom (max)
def is_pareto_efficient(costs, maximize_last=True):
    """
    Finds Pareto-efficient points given costs matrix where columns 0..-2 are minimized,
    and column -1 is maximized.
    Point j is dominated if there exists point i such that i is <= j on all min columns,
    >= j on max column, and strictly better on at least one.
    """
    n_points = costs.shape[0]
    is_efficient = np.ones(n_points, dtype=bool)
    for i in range(n_points):
        for j in range(n_points):
            if i == j or not is_efficient[j]:
                continue
            # Does point i dominate point j?
            if maximize_last:
                better_or_equal = np.all(costs[i, :-1] <= costs[j, :-1]) and (costs[i, -1] >= costs[j, -1])
                strictly_better = np.any(costs[i, :-1] < costs[j, :-1]) or (costs[i, -1] > costs[j, -1])
            else:
                better_or_equal = np.all(costs[i] <= costs[j])
                strictly_better = np.any(costs[i] < costs[j])
            if better_or_equal and strictly_better:
                is_efficient[j] = False
    return is_efficient

costs_nom = np.column_stack([
    df_nom_feas["Mass_kg"].values,
    df_nom_feas["Delta_P_nom_kPa"].values,
    df_nom_feas["T_clad_nom_K"].values,
    df_nom_feas["Lifetime_nom_yr"].values
])

pareto_mask_nom = is_pareto_efficient(costs_nom, maximize_last=True)
df_pareto_nom = df_nom_feas[pareto_mask_nom].copy()
df_pareto_nom.to_csv(os.path.join(OPT_DIR, "nominal_pareto_front.csv"), index=False)
print(f"Saved nominal_pareto_front.csv ({len(df_pareto_nom)} candidates on nominal front).")

# -----------------------------------------------------------------------------
# 3. IDENTIFY ROBUST PARETO FRONTIER
# -----------------------------------------------------------------------------
# Feasible robustly: Feasibility_Robust == "ROBUST PASS"
df_rob_feas = df_all[df_all["Feasibility_Robust"] == "ROBUST PASS"].copy()

costs_rob = np.column_stack([
    df_rob_feas["Mass_kg"].values,
    df_rob_feas["Delta_P_pert_kPa"].values,
    df_rob_feas["T_clad_pert_K"].values,
    df_rob_feas["Lifetime_robust_yr"].values
])

pareto_mask_rob = is_pareto_efficient(costs_rob, maximize_last=True)
df_pareto_rob = df_rob_feas[pareto_mask_rob].copy()
df_pareto_rob.to_csv(os.path.join(OPT_DIR, "robust_pareto_front.csv"), index=False)
print(f"Saved robust_pareto_front.csv ({len(df_pareto_rob)} candidates on robust front).")

print("\n=== NOMINAL PARETO FRONT ===")
print(df_pareto_nom[["Candidate_ID", "Mass_kg", "Delta_P_nom_kPa", "T_clad_nom_K", "Lifetime_nom_yr"]].to_string(index=False))

print("\n=== ROBUST PARETO FRONT ===")
print(df_pareto_rob[["Candidate_ID", "Mass_kg", "Delta_P_pert_kPa", "T_clad_pert_K", "Lifetime_robust_yr"]].to_string(index=False))
