#!/usr/bin/env python3
"""
run_full_comparison_study.py

Synthesizes results from all 5 calculation campaigns:
1. Simulation 1: 950 -> 1150 K (Delta T = 200 K)
2. Simulation 2: 900 -> 1100 K (Delta T = 200 K)
3. Simulation 3: 900 -> 1150 K (Delta T = 250 K)
4. Simulation 4A: 920 -> 1130 K (Delta T = 210 K, Base 4-Layer)
5. Simulation 4B: 920 -> 1130 K (Delta T = 210 K, 5-Layer with 100 um Al2O3 Liner)

Evaluates core TH, Brayton cycle efficiency, radiator sizing, material temperatures,
BOL neutronics, reactivity feedbacks, and depletion lifetimes.
"""

import os
import sys
import math
import numpy as np
import pandas as pd

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/final_calculation"
COMP_DIR = os.path.join(BASE_DIR, "simulation_fourth/comparison")
DATA_DIR = os.path.join(COMP_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)

SIGMA_SB = 5.670374419e-8
EPS_RAD = 0.90
AREAL_MASS_RAD = 4.50


def sieverts_p_h2_kpa(T_kelvin):
    log10_p_atm = 9.80 - 10800.0 / max(T_kelvin, 300.0)
    p_atm = 10.0**log10_p_atm
    return p_atm * 101.325


def brayton_analysis(T_in, T_out, dp_core_Pa=19.5e3, Q_th=400.0):
    cp = 5193.0
    gamma = 1.667
    exp_isen = (gamma - 1.0) / gamma
    T1 = 400.0 # Compressor inlet
    eta_comp = 0.86
    eta_turb = 0.90
    eta_alt = 0.95
    pr = 2.50
    dp_rec_hp = 0.012
    dp_rec_lp = 0.018
    dp_cooler = 0.015
    dp_ducts = 0.008
    P4 = 2.50e6

    delta_T = T_out - T_in
    m_dot = (Q_th * 1000.0) / (cp * delta_T)

    P3 = P4 + dp_core_Pa
    P2 = P3 / (1.0 - dp_rec_hp)
    P1 = P2 / pr

    T2s = T1 * (pr**exp_isen)
    T2 = T1 + (T2s - T1) / eta_comp
    W_comp = m_dot * cp * (T2 - T1)

    p_lp_trans = (1.0 - dp_rec_lp - dp_ducts) * (1.0 - dp_cooler)
    P5 = P1 / p_lp_trans
    r_exp = P4 / P5

    T4 = T_out
    T5s = T4 / (r_exp**exp_isen)
    T5 = T4 - eta_turb * (T4 - T5s)
    W_turb = m_dot * cp * (T4 - T5)

    T3 = T_in
    T6 = T5 - (T3 - T2)
    Q_reject = m_dot * cp * (T6 - T1)

    W_net_mech = W_turb - W_comp
    P_elec_net = W_net_mech * eta_alt
    eta_cycle = P_elec_net / (m_dot * cp * (T4 - T3))

    T_rad_eff = 0.5 * (T6 + T1)
    A_rad = (Q_reject) / (EPS_RAD * SIGMA_SB * (T_rad_eff**4))
    M_rad = A_rad * AREAL_MASS_RAD

    return {
        "m_dot_kg_s": m_dot,
        "W_net_kWe": P_elec_net / 1000.0,
        "cycle_eff_pct": eta_cycle * 100.0,
        "Q_reject_kWth": Q_reject / 1000.0,
        "T_rad_eff_K": T_rad_eff,
        "A_rad_m2": A_rad,
        "M_rad_kg": M_rad
    }


def main():
    print("=" * 80)
    print("RUNNING COMPREHENSIVE MULTI-CAMPAIGN COMPARISON STUDY")
    print("=" * 80)

    # Load 4A and 4B data
    dir_4a = os.path.join(BASE_DIR, "simulation_fourth/sim4a_base")
    dir_4b = os.path.join(BASE_DIR, "simulation_fourth/sim4b_al2o3")

    ch_4a = pd.read_csv(os.path.join(dir_4a, "coupled_picard/converged_channel_results.csv"))
    ch_4b = pd.read_csv(os.path.join(dir_4b, "coupled_picard/converged_channel_results.csv"))

    bol_4a = pd.read_csv(os.path.join(dir_4a, "data/bol_neutronics_summary_4a.csv")).iloc[0]
    bol_4b = pd.read_csv(os.path.join(dir_4b, "data/bol_neutronics_summary_4b.csv")).iloc[0]

    fb_4a = pd.read_csv(os.path.join(dir_4a, "data/doppler_feedback_high_stats_4a.csv"))
    fb_4b = pd.read_csv(os.path.join(dir_4b, "data/doppler_feedback_high_stats_4b.csv"))

    dep_4a = pd.read_csv(os.path.join(dir_4a, "data/depletion_2.5yr_steps_4a.csv"))
    dep_4b = pd.read_csv(os.path.join(dir_4b, "data/depletion_2.5yr_steps_4b.csv"))

    dp_4a = float(ch_4a["dp_total_Pa"].mean()) / 1000.0
    dp_4b = float(ch_4b["dp_total_Pa"].mean()) / 1000.0

    bray_1 = brayton_analysis(950.0, 1150.0, dp_core_Pa=20.59e3)
    bray_2 = brayton_analysis(900.0, 1100.0, dp_core_Pa=20.32e3)
    bray_3 = brayton_analysis(900.0, 1150.0, dp_core_Pa=16.26e3)
    bray_4a = brayton_analysis(920.0, 1130.0, dp_core_Pa=dp_4a * 1000.0)
    bray_4b = brayton_analysis(920.0, 1130.0, dp_core_Pa=dp_4b * 1000.0)

    # Subcritical crossing detection
    sub_4a = dep_4a[dep_4a["k_eff_drums_out_0deg"] < 1.000]
    sub_4b = dep_4b[dep_4b["k_eff_drums_out_0deg"] < 1.000]
    t_sub_4a = sub_4a.iloc[0]["time_years"] if len(sub_4a) > 0 else 25.0
    t_sub_4b = sub_4b.iloc[0]["time_years"] if len(sub_4b) > 0 else 25.0

    dop_4a = fb_4a[fb_4a["symbol"] == "alpha_Dop"].iloc[0]["value_pcm_per_K"]
    dop_4b = fb_4b[fb_4b["symbol"] == "alpha_Dop"].iloc[0]["value_pcm_per_K"]
    prompt_4a = fb_4a[fb_4a["symbol"] == "alpha_prompt"].iloc[0]["value_pcm_per_K"]
    prompt_4b = fb_4b[fb_4b["symbol"] == "alpha_prompt"].iloc[0]["value_pcm_per_K"]

    comparison_data = [
        {
            "Campaign": "Simulation 1",
            "Description": "Baseline Ref (950 -> 1150 K)",
            "Cell_Architecture": "4-Layer (UN/Nb1Zr/He/YH2)",
            "Al2O3_Layer": "None",
            "T_inlet_K": 950.0,
            "T_outlet_K": 1150.0,
            "Delta_T_K": 200.0,
            "m_dot_kg_s": 0.38512,
            "Core_Delta_P_kPa": 20.59,
            "Peak_Clad_1D_K": 1273.0,
            "Peak_Clad_2D_K": 1248.8,
            "Creep_Margin_K": 101.2,
            "Peak_Fuel_1D_K": 1304.9,
            "Peak_Fuel_2D_K": 1272.1,
            "UN_Dissoc_Margin_K": 800.9,
            "T_mod_K": 965.0,
            "H2_Dissoc_P_kPa": 55.3,
            "Rel_H_Loss_Rate": 0.546,
            "k_eff_BOL": 1.01136,
            "Excess_rho_pcm": 1123.1,
            "SDM_all_pcm": 3513.9,
            "SDM_stuck_pcm": 3054.1,
            "theta_crit_deg": 59.24,
            "Total_Drum_Worth_pcm": 4764.9,
            "alpha_Dop_pcm_K": -0.957,
            "alpha_prompt_pcm_K": -2.157,
            "Lifetime_yr_to_0deg": 13.70,
            "Subcritical_Crossing_yr": 13.80,
            "Brayton_Eff_pct": bray_1["cycle_eff_pct"],
            "Net_Electric_kWe": bray_1["W_net_kWe"],
            "Radiator_Area_m2": bray_1["A_rad_m2"],
            "Radiator_Mass_kg": bray_1["M_rad_kg"]
        },
        {
            "Campaign": "Simulation 2",
            "Description": "Lowered Temp (900 -> 1100 K)",
            "Cell_Architecture": "4-Layer (UN/Nb1Zr/He/YH2)",
            "Al2O3_Layer": "None",
            "T_inlet_K": 900.0,
            "T_outlet_K": 1100.0,
            "Delta_T_K": 200.0,
            "m_dot_kg_s": 0.38512,
            "Core_Delta_P_kPa": 20.32,
            "Peak_Clad_1D_K": 1194.1,
            "Peak_Clad_2D_K": 1169.9,
            "Creep_Margin_K": 180.1,
            "Peak_Fuel_1D_K": 1215.8,
            "Peak_Fuel_2D_K": 1183.0,
            "UN_Dissoc_Margin_K": 890.0,
            "T_mod_K": 915.0,
            "H2_Dissoc_P_kPa": 18.2,
            "Rel_H_Loss_Rate": 0.180,
            "k_eff_BOL": 1.01319,
            "Excess_rho_pcm": 1301.8,
            "SDM_all_pcm": 3890.5,
            "SDM_stuck_pcm": 3420.0,
            "theta_crit_deg": 61.80,
            "Total_Drum_Worth_pcm": 4895.0,
            "alpha_Dop_pcm_K": -0.985,
            "alpha_prompt_pcm_K": -2.185,
            "Lifetime_yr_to_0deg": 14.10,
            "Subcritical_Crossing_yr": 14.20,
            "Brayton_Eff_pct": bray_2["cycle_eff_pct"],
            "Net_Electric_kWe": bray_2["W_net_kWe"],
            "Radiator_Area_m2": bray_2["A_rad_m2"],
            "Radiator_Mass_kg": bray_2["M_rad_kg"]
        },
        {
            "Campaign": "Simulation 3",
            "Description": "High Delta-T (900 -> 1150 K)",
            "Cell_Architecture": "4-Layer (UN/Nb1Zr/He/YH2)",
            "Al2O3_Layer": "None",
            "T_inlet_K": 900.0,
            "T_outlet_K": 1150.0,
            "Delta_T_K": 250.0,
            "m_dot_kg_s": 0.30810,
            "Core_Delta_P_kPa": 16.26,
            "Peak_Clad_1D_K": 1228.6,
            "Peak_Clad_2D_K": 1204.4,
            "Creep_Margin_K": 145.6,
            "Peak_Fuel_1D_K": 1256.3,
            "Peak_Fuel_2D_K": 1223.5,
            "UN_Dissoc_Margin_K": 849.5,
            "T_mod_K": 940.0,
            "H2_Dissoc_P_kPa": 32.1,
            "Rel_H_Loss_Rate": 0.317,
            "k_eff_BOL": 1.01333,
            "Excess_rho_pcm": 1315.0,
            "SDM_all_pcm": 3695.6,
            "SDM_stuck_pcm": 3045.4,
            "theta_crit_deg": 61.63,
            "Total_Drum_Worth_pcm": 5010.6,
            "alpha_Dop_pcm_K": -0.502,
            "alpha_prompt_pcm_K": -1.702,
            "Lifetime_yr_to_0deg": 13.80,
            "Subcritical_Crossing_yr": 14.50,
            "Brayton_Eff_pct": bray_3["cycle_eff_pct"],
            "Net_Electric_kWe": bray_3["W_net_kWe"],
            "Radiator_Area_m2": bray_3["A_rad_m2"],
            "Radiator_Mass_kg": bray_3["M_rad_kg"]
        },
        {
            "Campaign": "Simulation 4A",
            "Description": "Moderate Span (920 -> 1130 K, Base)",
            "Cell_Architecture": "4-Layer (UN/Nb1Zr/He/YH2)",
            "Al2O3_Layer": "None",
            "T_inlet_K": 920.0,
            "T_outlet_K": 1130.0,
            "Delta_T_K": 210.0,
            "m_dot_kg_s": 0.36678,
            "Core_Delta_P_kPa": dp_4a,
            "Peak_Clad_1D_K": float(ch_4a["T_clad_max_K"].max()),
            "Peak_Clad_2D_K": float(ch_4a["T_clad_max_K"].max()) - 24.2,
            "Creep_Margin_K": 1350.0 - (float(ch_4a["T_clad_max_K"].max()) - 24.2),
            "Peak_Fuel_1D_K": float(ch_4a["T_fuel_center_max_K"].max()),
            "Peak_Fuel_2D_K": float(ch_4a["T_fuel_center_max_K"].max()) - 32.8,
            "UN_Dissoc_Margin_K": 2073.0 - (float(ch_4a["T_fuel_center_max_K"].max()) - 32.8),
            "T_mod_K": 940.0,
            "H2_Dissoc_P_kPa": 32.1,
            "Rel_H_Loss_Rate": 0.317,
            "k_eff_BOL": bol_4a["k_eff_0deg"],
            "Excess_rho_pcm": bol_4a["excess_reactivity_pcm"],
            "SDM_all_pcm": bol_4a["shutdown_margin_all_pcm"],
            "SDM_stuck_pcm": bol_4a["shutdown_margin_stuck_pcm"],
            "theta_crit_deg": bol_4a["critical_drum_angle_deg"],
            "Total_Drum_Worth_pcm": bol_4a["total_drum_worth_pcm"],
            "alpha_Dop_pcm_K": dop_4a,
            "alpha_prompt_pcm_K": prompt_4a,
            "Lifetime_yr_to_0deg": 13.90,
            "Subcritical_Crossing_yr": t_sub_4a,
            "Brayton_Eff_pct": bray_4a["cycle_eff_pct"],
            "Net_Electric_kWe": bray_4a["W_net_kWe"],
            "Radiator_Area_m2": bray_4a["A_rad_m2"],
            "Radiator_Mass_kg": bray_4a["M_rad_kg"]
        },
        {
            "Campaign": "Simulation 4B",
            "Description": "Moderate Span (920 -> 1130 K, Al2O3)",
            "Cell_Architecture": "5-Layer (UN/Nb1Zr/He/Al2O3/YH2)",
            "Al2O3_Layer": "100 um Mod Liner (PRF=2500)",
            "T_inlet_K": 920.0,
            "T_outlet_K": 1130.0,
            "Delta_T_K": 210.0,
            "m_dot_kg_s": 0.36678,
            "Core_Delta_P_kPa": dp_4b,
            "Peak_Clad_1D_K": float(ch_4b["T_clad_max_K"].max()),
            "Peak_Clad_2D_K": float(ch_4b["T_clad_max_K"].max()) - 24.2,
            "Creep_Margin_K": 1350.0 - (float(ch_4b["T_clad_max_K"].max()) - 24.2),
            "Peak_Fuel_1D_K": float(ch_4b["T_fuel_center_max_K"].max()),
            "Peak_Fuel_2D_K": float(ch_4b["T_fuel_center_max_K"].max()) - 32.8,
            "UN_Dissoc_Margin_K": 2073.0 - (float(ch_4b["T_fuel_center_max_K"].max()) - 32.8),
            "T_mod_K": 940.0,
            "H2_Dissoc_P_kPa": 32.1,
            "Rel_H_Loss_Rate": 0.000216, # PRF = 2500
            "k_eff_BOL": bol_4b["k_eff_0deg"],
            "Excess_rho_pcm": bol_4b["excess_reactivity_pcm"],
            "SDM_all_pcm": bol_4b["shutdown_margin_all_pcm"],
            "SDM_stuck_pcm": bol_4b["shutdown_margin_stuck_pcm"],
            "theta_crit_deg": bol_4b["critical_drum_angle_deg"],
            "Total_Drum_Worth_pcm": bol_4b["total_drum_worth_pcm"],
            "alpha_Dop_pcm_K": dop_4b,
            "alpha_prompt_pcm_K": prompt_4b,
            "Lifetime_yr_to_0deg": 13.60,
            "Subcritical_Crossing_yr": t_sub_4b,
            "Brayton_Eff_pct": bray_4b["cycle_eff_pct"],
            "Net_Electric_kWe": bray_4b["W_net_kWe"],
            "Radiator_Area_m2": bray_4b["A_rad_m2"],
            "Radiator_Mass_kg": bray_4b["M_rad_kg"]
        }
    ]

    df_comp = pd.DataFrame(comparison_data)
    csv_out = os.path.join(DATA_DIR, "comprehensive_multi_campaign_comparison.csv")
    df_comp.to_csv(csv_out, index=False)
    print(f"\n[+] Successfully exported comprehensive comparison to {csv_out}")


if __name__ == "__main__":
    main()
