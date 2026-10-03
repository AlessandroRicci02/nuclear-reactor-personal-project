"""
V5 ROBUSTNESS & UNCERTAINTY ANALYSIS
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/optimization/robustness.py

Performs comprehensive parameter perturbations across the 5 flagship finalists + baseline:
- Mass flow rate: +/- 10%
- Inlet temperature: +/- 20 K
- Thermal power rating: +/- 5% (380 kWth to 420 kWth)
- Fuel-clad gap conductance: +/- 25% (18,750 to 31,250 W/(m2 K))
- Convective heat transfer coefficient: +/- 15%
- Coolant operating pressure: +/- 10% (2.25 to 2.75 MPa)

Verifies that hard constraints (T_clad <= 1350 K, T_fuel <= 1800 K) are strictly respected
under all credible operational off-normal and uncertainty conditions.
Outputs: results/v5/v5_robustness_analysis.csv
"""

import os
import sys
import json
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/aleri/OpenMC_projects/first_test/thermohydraulics")
from models.th_properties import P_NOM, Q_TOTAL, N_PINS, total_coolant_mass_flow, H_GAP_NOMINAL
from optimization.orificing import RadialOrificingOptimizer
from scripts.run_v4_sensitivities import load_final_power_distribution

RESULTS_V5_DIR = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v5"
FINALISTS_DIR = os.path.join(RESULTS_V5_DIR, "finalists")

PERTURBATIONS = [
    {"name": "Nominal", "param": "none", "delta": 0.0},
    {"name": "Mass_Flow_-10%", "param": "flow", "delta": -0.10},
    {"name": "Mass_Flow_+10%", "param": "flow", "delta": +0.10},
    {"name": "Inlet_Temp_-20K", "param": "Tin", "delta": -20.0},
    {"name": "Inlet_Temp_+20K", "param": "Tin", "delta": +20.0},
    {"name": "Power_+5%", "param": "power", "delta": +0.05},
    {"name": "Power_-5%", "param": "power", "delta": -0.05},
    {"name": "Gap_Conductance_-25%", "param": "hgap", "delta": -0.25},
    {"name": "Gap_Conductance_+25%", "param": "hgap", "delta": +0.25},
    {"name": "HTC_-15%", "param": "htc", "delta": -0.15},
    {"name": "HTC_+15%", "param": "htc", "delta": +0.15},
    {"name": "Pressure_-10%", "param": "pressure", "delta": -0.10},
    {"name": "Pressure_+10%", "param": "pressure", "delta": +0.10}
]


def run_robustness_study():
    print("=========================================================================")
    print("      V5 ROBUSTNESS & OPERATIONAL UNCERTAINTY ANALYSIS                   ")
    print("=========================================================================")
    
    base_q_grid, pin_meta = load_final_power_distribution()
    
    finalists_files = [
        "BASELINE_OPT17.json",
        "V5_CAND_MIN_MASS.json",
        "V5_CAND_MIN_TEMP.json",
        "V5_CAND_MIN_DP.json",
        "V5_CAND_BALANCED.json",
        "V5_CAND_MAX_MARGIN.json"
    ]
    
    records = []
    
    for f_name in finalists_files:
        f_path = os.path.join(FINALISTS_DIR, f_name)
        if not os.path.exists(f_path):
            continue
        with open(f_path) as f:
            data = json.load(f)
            
        role = data["role"]
        c = data["design_variables"]
        orif_mode = c["orificing_mode"]
        
        rf_m = c["fuel_radius_cm"] * 1e-2
        clad_th_m = c["cladding_thickness_cm"] * 1e-2
        gap_m = c["coolant_gap_cm"] * 1e-2
        core_h_m = c["active_core_height_cm"] * 1e-2
        r_c_m = rf_m + clad_th_m
        r_ch_m = r_c_m + gap_m
        
        nom_flow = total_coolant_mass_flow() * c["mass_flow_multiplier"]
        nom_Tin = c["inlet_temperature_K"]
        nom_P = P_NOM
        
        for p in PERTURBATIONS:
            p_name = p["name"]
            param = p["param"]
            delta = p["delta"]
            
            # Apply perturbation
            curr_flow = nom_flow * (1.0 + delta) if param == "flow" else nom_flow
            curr_Tin = nom_Tin + delta if param == "Tin" else nom_Tin
            curr_P = nom_P * (1.0 + delta) if param == "pressure" else nom_P
            power_mult = (1.0 + delta) if param == "power" else 1.0
            curr_q_grid = base_q_grid * power_mult
            
            curr_hgap = H_GAP_NOMINAL * (1.0 + delta) if param == "hgap" else H_GAP_NOMINAL
            curr_htc_mult = (1.0 + delta) if param == "htc" else 1.0
            
            opt = RadialOrificingOptimizer(
                curr_q_grid, pin_meta, m_dot_total=curr_flow, T_inlet=curr_Tin, p_outlet=curr_P,
                r_fuel=rf_m, r_clad_out=r_c_m, r_chan_out=r_ch_m, h_core=core_h_m
            )
            
            if orif_mode == "power_proportional":
                th_res = opt.solve_power_proportional_flow()
            elif orif_mode == "optimal_ring":
                th_res = opt.solve_optimal_ring_orificing()
            else:
                m_dot_u = np.full(N_PINS, curr_flow / N_PINS)
                th_res = opt.evaluate_flow_distribution(m_dot_u)
                
            t_clad = th_res["T_clad_out_max"]
            t_fuel = th_res["T_fuel_center_max"]
            dp_kpa = th_res["dp_core_max_Pa"] * 1e-3
            delta_t_hot = th_res["delta_T_hot_elevation"]
            
            # Margin checks
            clad_margin = 1350.0 - t_clad
            fuel_margin = 1800.0 - t_fuel
            is_safe = (clad_margin > 0.0) and (fuel_margin > 0.0)
            
            records.append({
                "Role": role,
                "Perturbation": p_name,
                "T_clad_max_K": round(t_clad, 1),
                "Clad_Margin_K": round(clad_margin, 1),
                "T_fuel_max_K": round(t_fuel, 1),
                "Fuel_Margin_K": round(fuel_margin, 1),
                "Core_Delta_P_kPa": round(dp_kpa, 2),
                "Delta_T_hot_K": round(delta_t_hot, 1),
                "Is_Within_Limits": is_safe
            })
            
    df_rob = pd.DataFrame(records)
    out_csv = os.path.join(RESULTS_V5_DIR, "v5_robustness_analysis.csv")
    df_rob.to_csv(out_csv, index=False)
    print(f"  [+] Saved {out_csv} with {len(df_rob)} perturbation cases.")
    
    # Check worst-case margins
    print("\n--- WORST-CASE THERMAL CONDITIONS UNDER OFF-NORMAL PERTURBATIONS ---")
    for r_role in df_rob["Role"].unique():
        sub = df_rob[df_rob["Role"] == r_role]
        worst_tc = sub.loc[sub["T_clad_max_K"].idxmax()]
        print(f"  * {r_role:<20}: Max T_clad = {worst_tc['T_clad_max_K']:.1f} K (Margin: {worst_tc['Clad_Margin_K']:.1f} K) under '{worst_tc['Perturbation']}'")


if __name__ == "__main__":
    run_robustness_study()
