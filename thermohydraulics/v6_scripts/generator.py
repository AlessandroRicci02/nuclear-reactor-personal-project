"""
V5 DESIGN-SPACE GENERATOR & LEVEL 1 SCREENING DOE
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/optimization/generator.py

Generates and evaluates a population of candidate designs using Latin Hypercube Sampling (LHS)
plus targeted engineering candidates across the parameterized design space.
Checks all 9 hard engineering constraints and produces:
- results/v5/v5_screening_population.csv
- results/v5/candidates/V5_CAND_*.json
"""

import os
import sys
import json
import time
import math
import numpy as np
import pandas as pd
from scipy.stats import qmc

sys.path.insert(0, "/home/aleri/OpenMC_projects/first_test/thermohydraulics")
from optimization.design_space import (
    BASELINE_DESIGN_VECTOR, DESIGN_VARIABLES_BOUNDS,
    ENGINEERING_CONSTRAINTS, OBJECTIVE_FUNCTIONS
)
from optimization.evaluator import evaluate_candidate
from scripts.run_v4_sensitivities import load_final_power_distribution

RESULTS_V5_DIR = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v5"
CANDIDATES_DIR = os.path.join(RESULTS_V5_DIR, "candidates")
os.makedirs(CANDIDATES_DIR, exist_ok=True)


def generate_candidate_population(num_lhs=80, seed=42):
    """
    Generates a population of candidate design vectors:
    - Opt-17 baseline (Cand 17-Clean)
    - Targeted trade-off candidates (min mass, min temp, min dp, balanced, high margin)
    - Latin Hypercube Sampled (LHS) exploratory candidates
    """
    candidates = []
    
    # 1. Baseline Candidate
    c_base = {
        "candidate_id": "V5_CAND_BASE",
        "category": "Baseline",
        "description": "Reconstructed Opt-17 baseline (Cand 17-Clean, un-orificed)",
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
    candidates.append(c_base)

    # 2. Baseline with optimal 9-ring orificing
    c_base_orif = {
        "candidate_id": "V5_CAND_BASE_ORIF",
        "category": "Baseline+Orificed",
        "description": "Baseline Opt-17 geometry with optimal 9-ring flow orificing",
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
        "orificing_mode": "optimal_ring"
    }
    candidates.append(c_base_orif)

    # 3. Targeted engineering variations
    targeted_configs = [
        # Low mass focused
        {
            "candidate_id": "V5_TARGET_LOW_MASS_1",
            "category": "Targeted_LowMass",
            "description": "Compact core height and slimmer reflector to minimize mass",
            "fuel_radius_cm": 0.575,
            "cladding_thickness_cm": 0.045,
            "coolant_gap_cm": 0.075,
            "pin_pitch_cm": 2.220,
            "active_core_height_cm": 40.50,
            "radial_reflector_thickness_cm": 8.80,
            "mass_flow_multiplier": 1.02,
            "inlet_temperature_K": 950.0,
            "inner_enrichment_wt": 19.30,
            "outer_enrichment_wt": 19.75,
            "orificing_mode": "optimal_ring"
        },
        {
            "candidate_id": "V5_TARGET_LOW_MASS_2",
            "category": "Targeted_LowMass",
            "description": "Aggressive lightweight design within criticality limit",
            "fuel_radius_cm": 0.565,
            "cladding_thickness_cm": 0.040,
            "coolant_gap_cm": 0.070,
            "pin_pitch_cm": 2.200,
            "active_core_height_cm": 39.50,
            "radial_reflector_thickness_cm": 8.50,
            "mass_flow_multiplier": 1.05,
            "inlet_temperature_K": 945.0,
            "inner_enrichment_wt": 19.50,
            "outer_enrichment_wt": 19.75,
            "orificing_mode": "optimal_ring"
        },
        # Low temperature / high thermal margin focused
        {
            "candidate_id": "V5_TARGET_LOW_TEMP_1",
            "category": "Targeted_LowTemp",
            "description": "Enhanced flow and power proportional orificing for lowest peak T_clad",
            "fuel_radius_cm": 0.585,
            "cladding_thickness_cm": 0.050,
            "coolant_gap_cm": 0.080,
            "pin_pitch_cm": 2.280,
            "active_core_height_cm": 43.50,
            "radial_reflector_thickness_cm": 9.80,
            "mass_flow_multiplier": 1.10,
            "inlet_temperature_K": 930.0,
            "inner_enrichment_wt": 19.20,
            "outer_enrichment_wt": 19.75,
            "orificing_mode": "power_proportional"
        },
        {
            "candidate_id": "V5_TARGET_LOW_TEMP_2",
            "category": "Targeted_LowTemp",
            "description": "Higher flow multiplier and optimal ring orificing",
            "fuel_radius_cm": 0.580,
            "cladding_thickness_cm": 0.045,
            "coolant_gap_cm": 0.082,
            "pin_pitch_cm": 2.270,
            "active_core_height_cm": 42.45,
            "radial_reflector_thickness_cm": 10.00,
            "mass_flow_multiplier": 1.12,
            "inlet_temperature_K": 935.0,
            "inner_enrichment_wt": 19.00,
            "outer_enrichment_wt": 19.75,
            "orificing_mode": "optimal_ring"
        },
        # Low pressure drop focused
        {
            "candidate_id": "V5_TARGET_LOW_DP_1",
            "category": "Targeted_LowDP",
            "description": "Enlarged coolant gap (0.85 mm) to slash friction pressure drop",
            "fuel_radius_cm": 0.575,
            "cladding_thickness_cm": 0.045,
            "coolant_gap_cm": 0.088,
            "pin_pitch_cm": 2.300,
            "active_core_height_cm": 42.00,
            "radial_reflector_thickness_cm": 9.80,
            "mass_flow_multiplier": 0.98,
            "inlet_temperature_K": 950.0,
            "inner_enrichment_wt": 19.30,
            "outer_enrichment_wt": 19.75,
            "orificing_mode": "optimal_ring"
        },
        {
            "candidate_id": "V5_TARGET_LOW_DP_2",
            "category": "Targeted_LowDP",
            "description": "Wide coolant gap (0.92 mm) with uniform flow for minimum pumping loss",
            "fuel_radius_cm": 0.570,
            "cladding_thickness_cm": 0.045,
            "coolant_gap_cm": 0.092,
            "pin_pitch_cm": 2.320,
            "active_core_height_cm": 41.50,
            "radial_reflector_thickness_cm": 9.80,
            "mass_flow_multiplier": 0.96,
            "inlet_temperature_K": 950.0,
            "inner_enrichment_wt": 19.30,
            "outer_enrichment_wt": 19.75,
            "orificing_mode": "uniform"
        },
        # High reactivity / safety margin focused
        {
            "candidate_id": "V5_TARGET_MAX_MARGIN_1",
            "category": "Targeted_MaxMargin",
            "description": "Thick radial reflector (11.0 cm) and robust clad for maximum safety margin",
            "fuel_radius_cm": 0.590,
            "cladding_thickness_cm": 0.055,
            "coolant_gap_cm": 0.078,
            "pin_pitch_cm": 2.300,
            "active_core_height_cm": 44.00,
            "radial_reflector_thickness_cm": 11.00,
            "mass_flow_multiplier": 1.05,
            "inlet_temperature_K": 945.0,
            "inner_enrichment_wt": 19.30,
            "outer_enrichment_wt": 19.75,
            "orificing_mode": "optimal_ring"
        },
        # Balanced trade-off
        {
            "candidate_id": "V5_TARGET_BALANCED_1",
            "category": "Targeted_Balanced",
            "description": "Harmonized parameters balancing mass, temperature, and pressure drop",
            "fuel_radius_cm": 0.578,
            "cladding_thickness_cm": 0.048,
            "coolant_gap_cm": 0.080,
            "pin_pitch_cm": 2.260,
            "active_core_height_cm": 41.80,
            "radial_reflector_thickness_cm": 9.50,
            "mass_flow_multiplier": 1.03,
            "inlet_temperature_K": 948.0,
            "inner_enrichment_wt": 19.25,
            "outer_enrichment_wt": 19.75,
            "orificing_mode": "optimal_ring"
        }
    ]
    candidates.extend(targeted_configs)

    # 4. Latin Hypercube Sampling (LHS) Exploration
    lhs_sampler = qmc.LatinHypercube(d=10, seed=seed)
    sample = lhs_sampler.random(n=num_lhs)
    
    # Continuous variable bounds: [min, max]
    var_ranges = [
        ("fuel_radius_cm", 0.530, 0.605),
        ("cladding_thickness_cm", 0.038, 0.060),
        ("coolant_gap_cm", 0.065, 0.095),
        ("pin_pitch_cm", 2.180, 2.360),
        ("active_core_height_cm", 39.50, 45.50),
        ("radial_reflector_thickness_cm", 8.20, 11.20),
        ("mass_flow_multiplier", 0.92, 1.12),
        ("inlet_temperature_K", 920.0, 970.0),
        ("inner_enrichment_wt", 18.00, 19.50),
        ("outer_enrichment_wt", 18.80, 19.75) # Strict non-proliferation upper bound
    ]
    
    lower_bounds = np.array([r[1] for r in var_ranges])
    upper_bounds = np.array([r[2] for r in var_ranges])
    scaled_sample = qmc.scale(sample, lower_bounds, upper_bounds)
    
    orif_modes = ["uniform", "optimal_ring", "power_proportional"]
    
    for i in range(num_lhs):
        vec = scaled_sample[i]
        c_id = f"V5_LHS_{i+1:03d}"
        
        mode = orif_modes[i % len(orif_modes)]
        
        cand = {
            "candidate_id": c_id,
            "category": "LHS_DOE",
            "description": f"Latin Hypercube Sample #{i+1}",
            "fuel_radius_cm": round(float(vec[0]), 4),
            "cladding_thickness_cm": round(float(vec[1]), 4),
            "coolant_gap_cm": round(float(vec[2]), 4),
            "pin_pitch_cm": round(float(vec[3]), 4),
            "active_core_height_cm": round(float(vec[4]), 2),
            "radial_reflector_thickness_cm": round(float(vec[5]), 2),
            "mass_flow_multiplier": round(float(vec[6]), 3),
            "inlet_temperature_K": round(float(vec[7]), 1),
            "inner_enrichment_wt": round(float(vec[8]), 2),
            "outer_enrichment_wt": min(round(float(vec[9]), 2), 19.75),
            "orificing_mode": mode
        }
        candidates.append(cand)

    return candidates


def run_screening_campaign(num_lhs=80):
    print("=========================================================================")
    print("      V5 CANDIDATE POPULATION GENERATION & LEVEL 1 SCREENING             ")
    print("=========================================================================")
    
    q_vol_grid, pin_meta = load_final_power_distribution()
    candidates = generate_candidate_population(num_lhs=num_lhs)
    print(f"  [+] Generated population of {len(candidates)} candidate reactor designs.")
    
    eval_results = []
    
    t0 = time.time()
    for idx, c in enumerate(candidates):
        cid = c["candidate_id"]
        res = evaluate_candidate(cid, c, q_vol_grid, pin_meta, level=1)
        
        # Compute violations
        violations = [k for k, v in res["constraints"].items() if not v.get("pass", True)]
        is_feas = len(violations) == 0
        res["is_feasible"] = is_feas
        res["violation_count"] = len(violations)
        res["violations"] = violations

        # Save individual candidate JSON
        c_json_path = os.path.join(CANDIDATES_DIR, f"{cid}.json")
        with open(c_json_path, "w") as f:
            json.dump(res, f, indent=2)
            
        # Compile flat record for screening CSV
        objs = res["objectives"]
        constr = res["constraints"]
        th = res["thermal_hydraulics"]
        neu = res["neutronics"]
        masses = res["mass"]
        
        eval_results.append({
            "Candidate_ID": cid,
            "Category": c.get("category", "DOE"),
            "Is_Feasible": is_feas,
            "Violation_Count": len(violations),
            "Violations": "; ".join(violations) if violations else "None",
            # Objectives
            "F1_Mass_Subassembly_kg": objs["F1_subassembly_mass_kg"],
            "F2_Peak_Clad_Temp_K": objs["F2_peak_cladding_temp_K"],
            "F3_Core_Pressure_Drop_kPa": objs["F3_core_pressure_drop_kPa"],
            "F4_Power_Peaking_Fq": objs["F4_total_peaking_Fq"],
            "F5_Hot_Channel_Elevation_K": objs["F5_hot_channel_elevation_K"],
            # Key Variables
            "Fuel_Radius_cm": c["fuel_radius_cm"],
            "Clad_Thickness_cm": c["cladding_thickness_cm"],
            "Coolant_Gap_cm": c["coolant_gap_cm"],
            "Pin_Pitch_cm": c["pin_pitch_cm"],
            "Active_Core_Height_cm": c["active_core_height_cm"],
            "Reflector_Thickness_cm": c["radial_reflector_thickness_cm"],
            "Mass_Flow_Mult": c["mass_flow_multiplier"],
            "Inlet_Temp_K": c["inlet_temperature_K"],
            "Inner_Enrichment_wt": c["inner_enrichment_wt"],
            "Outer_Enrichment_wt": c["outer_enrichment_wt"],
            "Orificing_Mode": c["orificing_mode"],
            # Constraints Values
            "T_fuel_max_K": th["T_fuel_center_max_K"],
            "k_eff": neu["k_eff"],
            "SDM_all_in_pcm": neu["sdm_all_in_pcm"],
            "SDM_stuck_pcm": neu["sdm_stuck_pcm"],
            "Total_Module_Mass_kg": masses["m_module_kg"],
            "Fuel_Mass_kg": masses["m_fuel_kg"],
            "U235_Mass_kg": masses["m_u235_kg"]
        })
        
        if (idx + 1) % 15 == 0 or idx == len(candidates) - 1:
            print(f"  ... Evaluated {idx+1}/{len(candidates)} candidates ({time.time() - t0:.1f} s elapsed)")

    df_pop = pd.DataFrame(eval_results)
    out_csv = os.path.join(RESULTS_V5_DIR, "v5_screening_population.csv")
    df_pop.to_csv(out_csv, index=False)
    print(f"\n  [+] Successfully saved {out_csv} with {len(df_pop)} candidates.")
    
    n_feasible = int(df_pop["Is_Feasible"].sum())
    n_infeasible = len(df_pop) - n_feasible
    print(f"  [+] Feasibility breakdown: {n_feasible} Feasible ({n_feasible/len(df_pop)*100:.1f}%), {n_infeasible} Infeasible ({n_infeasible/len(df_pop)*100:.1f}%)")
    
    return df_pop


if __name__ == "__main__":
    run_screening_campaign(num_lhs=80)
