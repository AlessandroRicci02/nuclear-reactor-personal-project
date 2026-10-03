"""
V5 LOCAL SENSITIVITY / DESIGN-SPACE MAPPING
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/optimization/local_sensitivity.py

Performs a systematic local design-of-experiments around Opt-17:
- Perturbs candidate variables [-10%, -5%, baseline, +5%, +10%]
- Evaluates: k_eff, reactivity, T_f,max, T_c,max, T_cool,max, dp, mass, F_q, SDM
- Calculates local derivatives and normalized sensitivities
- Outputs: results/v5/v5_local_sensitivities.csv
"""

import os
import sys
import copy
import pandas as pd
import numpy as np

sys.path.insert(0, "/home/aleri/OpenMC_projects/first_test/thermohydraulics")
from optimization.design_space import BASELINE_DESIGN_VECTOR
from optimization.evaluator import evaluate_candidate
from scripts.run_v4_sensitivities import load_final_power_distribution

RESULTS_V5_DIR = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v5"
os.makedirs(RESULTS_V5_DIR, exist_ok=True)


def run_local_sensitivity_study():
    print("=========================================================================")
    print("      V5 LOCAL SENSITIVITY & DESIGN-SPACE MAPPING CAMPAIGN               ")
    print("=========================================================================")
    
    q_vol_grid, pin_meta = load_final_power_distribution()
    
    # Baseline vector
    base_vars = {
        "fuel_radius_cm": 0.580,
        "cladding_thickness_cm": 0.050,
        "coolant_gap_cm": 0.075,
        "pin_pitch_cm": 2.270,
        "active_core_height_cm": 42.45,
        "radial_reflector_thickness_cm": 9.80,
        "mass_flow_multiplier": 1.00,
        "inlet_temperature_K": 950.0,
        "inner_enrichment_wt": 19.30,
        "outer_enrichment_wt": 19.75,
        "orificing_mode": "uniform"
    }
    
    base_eval = evaluate_candidate("V5_BASE", base_vars, q_vol_grid, pin_meta, level=1)
    print(f"  [+] Baseline Opt-17 Evaluated:")
    print(f"      Mass = {base_eval['objectives']['F1_subassembly_mass_kg']:.2f} kg")
    print(f"      T_c,max = {base_eval['objectives']['F2_peak_cladding_temp_K']:.1f} K")
    print(f"      T_f,max = {base_eval['thermal_hydraulics']['T_fuel_center_max_K']:.1f} K")
    print(f"      Delta P = {base_eval['objectives']['F3_core_pressure_drop_kPa']:.2f} kPa")
    print(f"      k_eff = {base_eval['neutronics']['k_eff']:.5f}")

    # Variables to map with perturbations [-10%, -5%, +5%, +10%]
    perturb_specs = [
        {"var": "fuel_radius_cm", "name": "Fuel Radius", "nominal": 0.580, "deltas": [-0.10, -0.05, 0.05, 0.10], "units": "cm"},
        {"var": "cladding_thickness_cm", "name": "Clad Thickness", "nominal": 0.050, "deltas": [-0.10, -0.05, 0.05, 0.10], "units": "cm"},
        {"var": "coolant_gap_cm", "name": "Coolant Gap", "nominal": 0.075, "deltas": [-0.10, -0.05, 0.05, 0.10], "units": "cm"},
        {"var": "pin_pitch_cm", "name": "Pin Pitch", "nominal": 2.270, "deltas": [-0.05, -0.025, 0.025, 0.05], "units": "cm"},
        {"var": "active_core_height_cm", "name": "Active Core Height", "nominal": 42.45, "deltas": [-0.10, -0.05, 0.05, 0.10], "units": "cm"},
        {"var": "radial_reflector_thickness_cm", "name": "Reflector Thickness", "nominal": 9.80, "deltas": [-0.10, -0.05, 0.05, 0.10], "units": "cm"},
        {"var": "mass_flow_multiplier", "name": "Mass Flow Multiplier", "nominal": 1.00, "deltas": [-0.10, -0.05, 0.05, 0.10], "units": "-"},
        {"var": "inlet_temperature_K", "name": "Inlet Temperature", "nominal": 950.0, "deltas": [-0.05, -0.025, 0.025, 0.05], "units": "K"},
        {"var": "outer_enrichment_wt", "name": "Outer Enrichment", "nominal": 19.75, "deltas": [-0.06, -0.03, -0.015], "units": "wt%"} # Cannot exceed 19.75%
    ]

    records = []
    
    # Store baseline record
    records.append({
        "Variable": "Baseline",
        "Perturbation_Percent": 0.0,
        "Variable_Value": 0.0,
        "Mass_kg": base_eval["objectives"]["F1_subassembly_mass_kg"],
        "Delta_Mass_kg": 0.0,
        "T_clad_max_K": base_eval["objectives"]["F2_peak_cladding_temp_K"],
        "Delta_T_clad_K": 0.0,
        "T_fuel_max_K": base_eval["thermal_hydraulics"]["T_fuel_center_max_K"],
        "Delta_T_fuel_K": 0.0,
        "Pressure_Drop_kPa": base_eval["objectives"]["F3_core_pressure_drop_kPa"],
        "Delta_Pressure_Drop_kPa": 0.0,
        "k_eff": base_eval["neutronics"]["k_eff"],
        "Delta_k_pcm": 0.0,
        "F_q": base_eval["objectives"]["F4_total_peaking_Fq"],
        "Delta_F_q": 0.0,
        "SDM_all_in_pcm": base_eval["neutronics"]["sdm_all_in_pcm"]
    })

    derivatives = []

    for spec in perturb_specs:
        v_name = spec["var"]
        nom = spec["nominal"]
        
        # Evaluate perturbations
        for d in spec["deltas"]:
            val = nom * (1.0 + d) if v_name != "outer_enrichment_wt" else nom + d * nom
            # Clamp enrichment to non-proliferation limit
            if v_name == "outer_enrichment_wt" and val > 19.75:
                val = 19.75
                
            test_vars = copy.deepcopy(base_vars)
            test_vars[v_name] = val
            
            cand_res = evaluate_candidate(f"SENS_{v_name}_{d:+.2f}", test_vars, q_vol_grid, pin_meta, level=1)
            
            d_mass = cand_res["objectives"]["F1_subassembly_mass_kg"] - base_eval["objectives"]["F1_subassembly_mass_kg"]
            d_tc = cand_res["objectives"]["F2_peak_cladding_temp_K"] - base_eval["objectives"]["F2_peak_cladding_temp_K"]
            d_tf = cand_res["thermal_hydraulics"]["T_fuel_center_max_K"] - base_eval["thermal_hydraulics"]["T_fuel_center_max_K"]
            d_dp = cand_res["objectives"]["F3_core_pressure_drop_kPa"] - base_eval["objectives"]["F3_core_pressure_drop_kPa"]
            d_k_pcm = (cand_res["neutronics"]["k_eff"] - base_eval["neutronics"]["k_eff"]) * 1e5
            d_fq = cand_res["objectives"]["F4_total_peaking_Fq"] - base_eval["objectives"]["F4_total_peaking_Fq"]
            
            records.append({
                "Variable": spec["name"],
                "Perturbation_Percent": round(d * 100, 1),
                "Variable_Value": round(val, 4),
                "Mass_kg": cand_res["objectives"]["F1_subassembly_mass_kg"],
                "Delta_Mass_kg": round(d_mass, 2),
                "T_clad_max_K": cand_res["objectives"]["F2_peak_cladding_temp_K"],
                "Delta_T_clad_K": round(d_tc, 2),
                "T_fuel_max_K": cand_res["thermal_hydraulics"]["T_fuel_center_max_K"],
                "Delta_T_fuel_K": round(d_tf, 2),
                "Pressure_Drop_kPa": cand_res["objectives"]["F3_core_pressure_drop_kPa"],
                "Delta_Pressure_Drop_kPa": round(d_dp, 3),
                "k_eff": cand_res["neutronics"]["k_eff"],
                "Delta_k_pcm": round(d_k_pcm, 1),
                "F_q": cand_res["objectives"]["F4_total_peaking_Fq"],
                "Delta_F_q": round(d_fq, 3),
                "SDM_all_in_pcm": cand_res["neutronics"]["sdm_all_in_pcm"]
            })

        # Calculate central difference derivatives around +/- 5% or available bounds
        d_pos = max(spec["deltas"])
        d_neg = min(spec["deltas"])
        val_pos = nom * (1.0 + d_pos)
        val_neg = nom * (1.0 + d_neg)
        
        vars_pos = copy.deepcopy(base_vars)
        vars_neg = copy.deepcopy(base_vars)
        vars_pos[v_name] = val_pos
        vars_neg[v_name] = val_neg
        
        res_pos = evaluate_candidate("DP", vars_pos, q_vol_grid, pin_meta, level=1)
        res_neg = evaluate_candidate("DN", vars_neg, q_vol_grid, pin_meta, level=1)
        
        dx = val_pos - val_neg
        
        d_mass_dx = (res_pos["objectives"]["F1_subassembly_mass_kg"] - res_neg["objectives"]["F1_subassembly_mass_kg"]) / dx
        d_tc_dx = (res_pos["objectives"]["F2_peak_cladding_temp_K"] - res_neg["objectives"]["F2_peak_cladding_temp_K"]) / dx
        d_tf_dx = (res_pos["thermal_hydraulics"]["T_fuel_center_max_K"] - res_neg["thermal_hydraulics"]["T_fuel_center_max_K"]) / dx
        d_dp_dx = (res_pos["objectives"]["F3_core_pressure_drop_kPa"] - res_neg["objectives"]["F3_core_pressure_drop_kPa"]) / dx
        d_k_dx = (res_pos["neutronics"]["k_eff"] - res_neg["neutronics"]["k_eff"]) * 1e5 / dx # pcm / unit
        d_fq_dx = (res_pos["objectives"]["F4_total_peaking_Fq"] - res_neg["objectives"]["F4_total_peaking_Fq"]) / dx

        # Normalized sensitivities: (dY/Y) / (dX/X)
        norm_mass = (d_mass_dx * nom) / base_eval["objectives"]["F1_subassembly_mass_kg"]
        norm_tc = (d_tc_dx * nom) / base_eval["objectives"]["F2_peak_cladding_temp_K"]
        norm_dp = (d_dp_dx * nom) / base_eval["objectives"]["F3_core_pressure_drop_kPa"]
        norm_k = (d_k_dx * nom) / (base_eval["neutronics"]["k_eff"] * 1e5)
        
        derivatives.append({
            "Variable": spec["name"],
            "Units": spec["units"],
            "d(Mass)/dx": round(d_mass_dx, 2),
            "d(T_clad)/dx": round(d_tc_dx, 2),
            "d(T_fuel)/dx": round(d_tf_dx, 2),
            "d(Delta_P)/dx": round(d_dp_dx, 3),
            "d(k)/dx [pcm/unit]": round(d_k_dx, 1),
            "d(F_q)/dx": round(d_fq_dx, 3),
            "Norm_Sens_Mass": round(norm_mass, 3),
            "Norm_Sens_T_clad": round(norm_tc, 3),
            "Norm_Sens_Delta_P": round(norm_dp, 3),
            "Norm_Sens_k_eff": round(norm_k, 3)
        })

    df_records = pd.DataFrame(records)
    out_csv = os.path.join(RESULTS_V5_DIR, "v5_local_sensitivities.csv")
    df_records.to_csv(out_csv, index=False)
    print(f"  [+] Saved {out_csv} with {len(df_records)} perturbation points.")

    df_derivs = pd.DataFrame(derivatives)
    out_derivs_csv = os.path.join(RESULTS_V5_DIR, "v5_derivatives_summary.csv")
    df_derivs.to_csv(out_derivs_csv, index=False)
    print(f"  [+] Saved {out_derivs_csv} with normalized derivatives.")
    
    print("\n--- KEY CONTROLLING VARIABLES & NORMALIZED SENSITIVITIES ---")
    for _, r in df_derivs.iterrows():
        print(f"  * {r['Variable']:<22}: Norm Sens [Mass: {r['Norm_Sens_Mass']:+.2f} | T_clad: {r['Norm_Sens_T_clad']:+.2f} | Delta_P: {r['Norm_Sens_Delta_P']:+.2f} | k_eff: {r['Norm_Sens_k_eff']:+.3f}]")

if __name__ == "__main__":
    run_local_sensitivity_study()
