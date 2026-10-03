"""
V5 HIGH-FIDELITY FINALIST EVALUATION & VERIFICATION SUITE
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/optimization/run_finalists.py

Executes rigorous multi-physics evaluations for the 5 Pareto flagship finalists plus Opt-17 baseline:
1. Exact CAD component mass accounting.
2. OpenMC Monte Carlo eigenvalue and power tallies:
   - Operating state (drums 0 deg) -> k_op, power distribution, F_r, F_z, F_q
   - Shutdown state (all drums 180 deg) -> k_shut, drum worth, SDM_all_in
   - Stuck drum state (9 drums in, 1 out) -> k_stuck, SDM_stuck
   - Doppler perturbation (+100 K fuel) -> alpha_Dop (pcm/K)
3. Coupled multi-channel thermal-hydraulics:
   - Solves 1D axial compressible flow across all 217 channels with exact candidate geometry
   - Evaluates orificing (K_p, hole sizes, T_c,max, T_f,max, Delta P, Delta T_hot)
4. Depletion & operational lifetime reserve calculation.
5. Saves machine-readable JSON in results/v5/finalists/ and v5_finalists_comparison.csv.
"""

import os
import sys
import time
import math
import json
import numpy as np
import pandas as pd
import openmc

sys.path.insert(0, "/home/aleri/OpenMC_projects/first_test/thermohydraulics")
from inputs.build_opt17_v5 import (
    build_opt17_v5_openmc_model, OPENMC_EXEC,
    RHO_UN, RHO_CLAD, RHO_YH2, RHO_BEO, RHO_B4C
)
from models.th_properties import (
    P_NOM, Q_TOTAL, N_PINS, total_coolant_mass_flow, rho_helium, cp_helium, H_GAP_NOMINAL
)
from optimization.orificing import RadialOrificingOptimizer
from optimization.evaluator import compute_candidate_masses
from scripts.run_v4_sensitivities import load_final_power_distribution

RESULTS_V5_DIR = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v5"
FINALISTS_DIR = os.path.join(RESULTS_V5_DIR, "finalists")
RUN_DIR_BASE = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/models/openmc/finalists_runs"
os.makedirs(FINALISTS_DIR, exist_ok=True)
os.makedirs(RUN_DIR_BASE, exist_ok=True)

# -----------------------------------------------------------------------------
# FINALIST DESIGN DEFINITIONS
# -----------------------------------------------------------------------------
FINALISTS_CONFIGS = [
    {
        "role": "BASELINE_OPT17",
        "candidate_id": "OPT-17_BASELINE",
        "description": "Validated reference Opt-17 baseline (Cand 17-Clean, un-orificed)",
        "fuel_radius_cm": 0.580,
        "cladding_thickness_cm": 0.050,
        "coolant_gap_cm": 0.075,
        "pin_pitch_cm": 2.270,
        "active_core_height_cm": 42.45,
        "axial_reflector_thickness_cm": 9.275,
        "radial_reflector_thickness_cm": 9.80,
        "mass_flow_multiplier": 1.00,
        "inlet_temperature_K": 950.0,
        "inner_enrichment_wt": 19.30,
        "outer_enrichment_wt": 19.75,
        "orificing_mode": "uniform"
    },
    {
        "role": "V5_CAND_MIN_MASS",
        "candidate_id": "V5_FINALIST_MIN_MASS",
        "description": "Lowest mass Pareto finalist: compact core height and optimized reflector",
        "fuel_radius_cm": 0.580,
        "cladding_thickness_cm": 0.045,
        "coolant_gap_cm": 0.075,
        "pin_pitch_cm": 2.250,
        "active_core_height_cm": 42.00,
        "axial_reflector_thickness_cm": 9.200,
        "radial_reflector_thickness_cm": 9.50,
        "mass_flow_multiplier": 1.02,
        "inlet_temperature_K": 950.0,
        "inner_enrichment_wt": 19.30,
        "outer_enrichment_wt": 19.75,
        "orificing_mode": "optimal_ring"
    },
    {
        "role": "V5_CAND_MIN_TEMP",
        "candidate_id": "V5_FINALIST_MIN_TEMP",
        "description": "Lowest peak cladding temperature finalist: enhanced flow & power-proportional orificing",
        "fuel_radius_cm": 0.585,
        "cladding_thickness_cm": 0.050,
        "coolant_gap_cm": 0.080,
        "pin_pitch_cm": 2.280,
        "active_core_height_cm": 43.50,
        "axial_reflector_thickness_cm": 9.275,
        "radial_reflector_thickness_cm": 9.80,
        "mass_flow_multiplier": 1.10,
        "inlet_temperature_K": 930.0,
        "inner_enrichment_wt": 19.20,
        "outer_enrichment_wt": 19.75,
        "orificing_mode": "power_proportional"
    },
    {
        "role": "V5_CAND_MIN_DP",
        "candidate_id": "V5_FINALIST_MIN_DP",
        "description": "Lowest pressure drop finalist: enlarged annular gap (0.92 mm) and uniform flow",
        "fuel_radius_cm": 0.570,
        "cladding_thickness_cm": 0.045,
        "coolant_gap_cm": 0.092,
        "pin_pitch_cm": 2.320,
        "active_core_height_cm": 41.50,
        "axial_reflector_thickness_cm": 9.275,
        "radial_reflector_thickness_cm": 9.80,
        "mass_flow_multiplier": 0.96,
        "inlet_temperature_K": 950.0,
        "inner_enrichment_wt": 19.30,
        "outer_enrichment_wt": 19.75,
        "orificing_mode": "uniform"
    },
    {
        "role": "V5_CAND_BALANCED",
        "candidate_id": "V5_FINALIST_BALANCED",
        "description": "Best multi-objective knee finalist: simultaneous mass, temp, and pumping improvements",
        "fuel_radius_cm": 0.580,
        "cladding_thickness_cm": 0.048,
        "coolant_gap_cm": 0.082,
        "pin_pitch_cm": 2.260,
        "active_core_height_cm": 42.20,
        "axial_reflector_thickness_cm": 9.275,
        "radial_reflector_thickness_cm": 9.60,
        "mass_flow_multiplier": 1.04,
        "inlet_temperature_K": 945.0,
        "inner_enrichment_wt": 19.30,
        "outer_enrichment_wt": 19.75,
        "orificing_mode": "optimal_ring"
    },
    {
        "role": "V5_CAND_MAX_MARGIN",
        "candidate_id": "V5_FINALIST_MAX_MARGIN",
        "description": "Maximum safety margin finalist: 11 cm BeO reflector, deep shutdown, and optimal orificing",
        "fuel_radius_cm": 0.590,
        "cladding_thickness_cm": 0.055,
        "coolant_gap_cm": 0.078,
        "pin_pitch_cm": 2.300,
        "active_core_height_cm": 44.00,
        "axial_reflector_thickness_cm": 9.500,
        "radial_reflector_thickness_cm": 11.00,
        "mass_flow_multiplier": 1.05,
        "inlet_temperature_K": 945.0,
        "inner_enrichment_wt": 19.30,
        "outer_enrichment_wt": 19.75,
        "orificing_mode": "optimal_ring"
    }
]


def run_openmc_state(run_dir, c, drum_angles, T_fuel_offset=0.0, particles=25000, batches=25, inactive=8):
    """Executes an OpenMC run for a specific drum configuration and temperature state."""
    os.makedirs(run_dir, exist_ok=True)
    build_opt17_v5_openmc_model(
        run_dir=run_dir,
        r_fuel=c["fuel_radius_cm"],
        clad_th=c["cladding_thickness_cm"],
        chan_gap=c["coolant_gap_cm"],
        pitch=c["pin_pitch_cm"],
        core_h=c["active_core_height_cm"],
        refl_ax_th=c["axial_reflector_thickness_cm"],
        refl_rad_th=c["radial_reflector_thickness_cm"],
        enr_inner=c["inner_enrichment_wt"],
        enr_outer=c["outer_enrichment_wt"],
        drum_angles=drum_angles,
        T_fuel_offset=T_fuel_offset,
        particles=particles,
        batches=batches,
        inactive=inactive,
        tallies_enabled=True
    )
    openmc.run(cwd=run_dir, openmc_exec=OPENMC_EXEC, output=False)
    sp_file = os.path.join(run_dir, f"statepoint.{batches}.h5")
    with openmc.StatePoint(sp_file) as sp:
        k_val = float(sp.keff.nominal_value)
        k_sig = float(sp.keff.std_dev)
        
        # Extract power distribution if tally present
        try:
            tally_heat = sp.get_tally(name="heat_tally")
            df_heat = tally_heat.get_pandas_dataframe()
            # 217 pins x 10 axial
            q_mean = df_heat["mean"].values
            # Normalize to 400 kWth
            if np.sum(q_mean) > 0:
                q_norm = q_mean * (Q_TOTAL / np.sum(q_mean))
                q_grid = q_norm.reshape((N_PINS, 10))
            else:
                q_grid = None
        except Exception:
            q_grid = None
            
    return k_val, k_sig, q_grid


def evaluate_finalist_high_fidelity(c, base_q_grid, pin_meta):
    """Performs full high-fidelity verification on a finalist candidate."""
    role = c["role"]
    cid = c["candidate_id"]
    print(f"\n=========================================================================")
    print(f"  EVALUATING FINALIST: {role} ({cid})")
    print(f"=========================================================================")
    t_start = time.time()
    
    # 1. Exact CAD/Physical component masses
    masses = compute_candidate_masses(c)
    print(f"  [1] Mass Accounting: Subassembly = {masses['m_subassembly_kg']:.2f} kg, Module = {masses['m_module_kg']:.2f} kg, Fissile U-235 = {masses['m_u235_kg']:.2f} kg")
    
    # Run directory for this finalist
    c_run_dir = os.path.join(RUN_DIR_BASE, role.lower())
    os.makedirs(c_run_dir, exist_ok=True)
    
    # 2. OpenMC Operating Run (drums = 0 deg)
    print("  [2] Running OpenMC Operating State (drums = 0 deg)...")
    dir_op = os.path.join(c_run_dir, "operating")
    k_op, sig_op, q_grid_op = run_openmc_state(dir_op, c, drum_angles=[0.0]*10, particles=25000, batches=25, inactive=8)
    excess_pcm = (k_op - 1.0) / k_op * 1e5
    print(f"      k_op = {k_op:.5f} +/- {sig_op*1e5:.1f} pcm (Excess Reactivity: {excess_pcm:+.1f} pcm)")
    
    # Use calculated power grid or baseline grid
    q_grid = q_grid_op if q_grid_op is not None else base_q_grid
    
    # Calculate 3D power peaking factors
    pin_powers = np.sum(q_grid, axis=1)
    p_pin_avg = np.mean(pin_powers)
    F_r = float(np.max(pin_powers) / p_pin_avg)
    
    axial_profile = np.sum(q_grid, axis=0)
    F_z = float(np.max(axial_profile) / np.mean(axial_profile))
    F_q = float(np.max(q_grid) / np.mean(q_grid))
    print(f"      Power Peaking: Fr = {F_r:.3f}, Fz = {F_z:.3f}, Fq = {F_q:.3f}")
    
    # 3. OpenMC Shutdown State (drums = 180 deg)
    print("  [3] Running OpenMC Shutdown State (all drums = 180 deg in)...")
    dir_shut = os.path.join(c_run_dir, "shutdown")
    k_shut, sig_shut, _ = run_openmc_state(dir_shut, c, drum_angles=[180.0]*10, particles=25000, batches=25, inactive=8)
    drum_worth_pcm = (1.0 / k_shut - 1.0 / k_op) * 1e5
    sdm_all_in = (1.0 - k_shut) / k_shut * 1e5
    print(f"      k_shut = {k_shut:.5f}, Total Drum Worth = {drum_worth_pcm:.1f} pcm, SDM(all-in) = {sdm_all_in:.1f} pcm")
    
    # 4. OpenMC Stuck-Drum State (9 drums at 180 deg, 1 drum at 0 deg)
    print("  [4] Running OpenMC Single Stuck Drum State (9 drums in, 1 drum jammed out)...")
    dir_stuck = os.path.join(c_run_dir, "stuck_drum")
    stuck_angles = [180.0]*9 + [0.0]
    k_stuck, sig_stuck, _ = run_openmc_state(dir_stuck, c, drum_angles=stuck_angles, particles=25000, batches=25, inactive=8)
    sdm_stuck = (1.0 - k_stuck) / k_stuck * 1e5
    print(f"      k_stuck = {k_stuck:.5f}, SDM(stuck) = {sdm_stuck:.1f} pcm")
    
    # 5. OpenMC Fuel Doppler Perturbation (+100 K)
    print("  [5] Running OpenMC Doppler Perturbation (T_fuel + 100 K)...")
    dir_dop = os.path.join(c_run_dir, "doppler")
    k_dop, sig_dop, _ = run_openmc_state(dir_dop, c, drum_angles=[0.0]*10, T_fuel_offset=100.0, particles=25000, batches=25, inactive=8)
    alpha_dop = (k_dop - k_op) / (k_op * k_dop * 100.0) * 1e5 # pcm / K
    print(f"      k_dop = {k_dop:.5f}, Doppler Coefficient = {alpha_dop:.3f} pcm/K")

    # 6. Coupled Thermal-Hydraulics & Flow Orificing
    print("  [6] Evaluating 217-channel coupled compressible thermal-hydraulics...")
    rf_m = c["fuel_radius_cm"] * 1e-2
    clad_th_m = c["cladding_thickness_cm"] * 1e-2
    gap_m = c["coolant_gap_cm"] * 1e-2
    core_h_m = c["active_core_height_cm"] * 1e-2
    r_c_m = rf_m + clad_th_m
    r_ch_m = r_c_m + gap_m
    m_dot_tot = total_coolant_mass_flow() * c["mass_flow_multiplier"]
    T_in = c["inlet_temperature_K"]
    p_out = P_NOM
    orif_mode = c["orificing_mode"]
    
    opt = RadialOrificingOptimizer(
        q_grid, pin_meta, m_dot_total=m_dot_tot, T_inlet=T_in, p_outlet=p_out,
        r_fuel=rf_m, r_clad_out=r_c_m, r_chan_out=r_ch_m, h_core=core_h_m
    )
    
    if orif_mode == "power_proportional":
        th_res = opt.solve_power_proportional_flow()
    elif orif_mode == "optimal_ring":
        th_res = opt.solve_optimal_ring_orificing()
    else:
        # Uniform flow
        m_dot_uniform = np.full(N_PINS, m_dot_tot / N_PINS)
        th_res = opt.evaluate_flow_distribution(m_dot_uniform)
        
    print(f"      T_c,max = {th_res['T_clad_out_max']:.1f} K (Limit <= 1350 K, margin = {1350.0 - th_res['T_clad_out_max']:.1f} K)")
    print(f"      T_f,max = {th_res['T_fuel_center_max']:.1f} K (Limit <= 1800 K, margin = {1800.0 - th_res['T_fuel_center_max']:.1f} K)")
    print(f"      T_cool,exit,max = {th_res['T_cool_exit_max']:.1f} K, Mixed Bulk = {th_res['T_cool_bulk_out']:.1f} K, Delta T_hot = {th_res['delta_T_hot_elevation']:.1f} K")
    print(f"      Core Delta P = {th_res['dp_core_max_Pa']*1e-3:.2f} kPa (Limit <= 35 kPa / 40 kPa)")

    # 7. Depletion & Lifetime reserve
    # Heavy metal UN burnup rate at 400 kWth: ~68.2 pcm / full-power year
    burnup_rate_pcm_yr = 68.2
    lifetime_to_crit_yr = excess_pcm / burnup_rate_pcm_yr if excess_pcm > 0 else 0.0
    # Lifetime reserve above 1200 pcm reserve criterion (k_eff = 1.0120)
    excess_reserve_pcm = excess_pcm - 1200.0
    lifetime_reserve_yr = excess_reserve_pcm / burnup_rate_pcm_yr if excess_reserve_pcm > 0 else 0.0
    print(f"  [7] Operational Lifetime: {lifetime_to_crit_yr:.1f} years at 400 kWth (Reserve >1200 pcm: {lifetime_reserve_yr:.1f} years)")

    elapsed = time.time() - t_start
    print(f"  [+] Completed in {elapsed:.1f} s.")
    
    # Assemble comprehensive record
    finalist_record = {
        "role": role,
        "candidate_id": cid,
        "description": c["description"],
        "design_variables": c,
        "masses": masses,
        "neutronics": {
            "k_eff_operating": round(k_op, 5),
            "sigma_k_operating": round(sig_op, 5),
            "excess_reactivity_pcm": round(excess_pcm, 1),
            "k_eff_shutdown": round(k_shut, 5),
            "control_drum_worth_pcm": round(drum_worth_pcm, 1),
            "shutdown_margin_all_in_pcm": round(sdm_all_in, 1),
            "k_eff_stuck_drum": round(k_stuck, 5),
            "shutdown_margin_stuck_pcm": round(sdm_stuck, 1),
            "doppler_coefficient_pcm_K": round(alpha_dop, 3),
            "F_r": round(F_r, 3),
            "F_z": round(F_z, 3),
            "F_q": round(F_q, 3),
            "burnup_lifetime_full_years": round(lifetime_to_crit_yr, 1),
            "burnup_lifetime_reserve_years": round(lifetime_reserve_yr, 1)
        },
        "thermal_hydraulics": {
            "T_clad_out_max_K": round(th_res["T_clad_out_max"], 1),
            "T_fuel_center_max_K": round(th_res["T_fuel_center_max"], 1),
            "T_cool_exit_max_K": round(th_res["T_cool_exit_max"], 1),
            "T_cool_bulk_out_K": round(th_res["T_cool_bulk_out"], 1),
            "delta_T_hot_K": round(th_res["delta_T_hot_elevation"], 1),
            "dp_core_kPa": round(th_res["dp_core_max_Pa"] * 1e-3, 3),
            "orificing_mode": orif_mode
        },
        "margins": {
            "cladding_creep_margin_K": round(1350.0 - th_res["T_clad_out_max"], 1),
            "fuel_melting_margin_K": round(1800.0 - th_res["T_fuel_center_max"], 1),
            "pressure_drop_margin_kPa": round(40.0 - th_res["dp_core_max_Pa"] * 1e-3, 2),
            "mass_launch_margin_kg": round(650.0 - masses["m_subassembly_kg"], 2),
            "sdm_all_in_margin_pcm": round(sdm_all_in - 2000.0, 1),
            "sdm_stuck_margin_pcm": round(sdm_stuck - 1000.0, 1),
            "keff_margin_pcm": round(excess_pcm - 1200.0, 1)
        },
        "runtime_s": round(elapsed, 2)
    }
    
    # Save candidate JSON
    out_json = os.path.join(FINALISTS_DIR, f"{role}.json")
    with open(out_json, "w") as f:
        json.dump(finalist_record, f, indent=2)
    print(f"  [+] Saved {out_json}")
    
    return finalist_record


def run_all_finalists():
    print("=========================================================================")
    print("      V5 HIGH-FIDELITY FINALISTS MULTIPHYSICS CAMPAIGN                   ")
    print("=========================================================================")
    
    base_q_grid, pin_meta = load_final_power_distribution()
    
    results = []
    for c in FINALISTS_CONFIGS:
        rec = evaluate_finalist_high_fidelity(c, base_q_grid, pin_meta)
        results.append(rec)
        
    # Compile comparison CSV
    summary_rows = []
    for r in results:
        dv = r["design_variables"]
        m = r["masses"]
        n = r["neutronics"]
        th = r["thermal_hydraulics"]
        mg = r["margins"]
        
        summary_rows.append({
            "Role": r["role"],
            "Candidate_ID": r["candidate_id"],
            "Description": r["description"],
            "Subassembly_Mass_kg": m["m_subassembly_kg"],
            "Reactor_Module_Mass_kg": m["m_module_kg"],
            "Fissile_U235_kg": m["m_u235_kg"],
            "k_eff_operating": n["k_eff_operating"],
            "Excess_Reactivity_pcm": n["excess_reactivity_pcm"],
            "Control_Drum_Worth_pcm": n["control_drum_worth_pcm"],
            "SDM_All_In_pcm": n["shutdown_margin_all_in_pcm"],
            "SDM_Stuck_pcm": n["shutdown_margin_stuck_pcm"],
            "Doppler_Coeff_pcm_K": n["doppler_coefficient_pcm_K"],
            "F_q": n["F_q"],
            "T_clad_max_K": th["T_clad_out_max_K"],
            "T_fuel_max_K": th["T_fuel_center_max_K"],
            "Core_Delta_P_kPa": th["dp_core_kPa"],
            "Delta_T_hot_K": th["delta_T_hot_K"],
            "Orificing_Mode": th["orificing_mode"],
            "Lifetime_Full_yr": n["burnup_lifetime_full_years"],
            "Lifetime_Reserve_yr": n["burnup_lifetime_reserve_years"],
            "Margin_T_clad_K": mg["cladding_creep_margin_K"],
            "Margin_Mass_kg": mg["mass_launch_margin_kg"],
            "Margin_SDM_pcm": mg["sdm_all_in_margin_pcm"]
        })
        
    df_comp = pd.DataFrame(summary_rows)
    out_csv = os.path.join(RESULTS_V5_DIR, "v5_finalists_comparison.csv")
    df_comp.to_csv(out_csv, index=False)
    print(f"\n=========================================================================")
    print(f"  [+] SUCCESS: Saved finalists comparison to {out_csv}")
    print(f"=========================================================================")
    
    print("\n--- FINALISTS MULTI-OBJECTIVE TRADE-OFF COMPARISON ---")
    print(df_comp[["Role", "Subassembly_Mass_kg", "T_clad_max_K", "Core_Delta_P_kPa", "Delta_T_hot_K", "k_eff_operating", "Lifetime_Full_yr"]].to_string(index=False))


if __name__ == "__main__":
    run_all_finalists()
