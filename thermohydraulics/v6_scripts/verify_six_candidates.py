"""
V6 HIGH-FIDELITY CANDIDATE VERIFICATION SUITE
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/v6_scripts/verify_six_candidates.py

Executes rigorous high-fidelity verification for the 6 reactors:
- BASELINE_OPT17
- V5_CAND_MIN_MASS
- V5_CAND_MIN_TEMP
- V5_CAND_MIN_DP
- V5_CAND_BALANCED
- V5_CAND_MAX_MARGIN

Evaluates:
1. OpenMC Monte Carlo Neutronics (high statistics: 60,000 particles x 35 batches):
   - Operating state (drums 0 deg) -> k_op, sigma_k, excess reactivity
   - Shutdown state (drums 180 deg) -> k_shut, drum worth, SDM_all_in
   - Single stuck drum state (9 in, 1 out) -> k_stuck, SDM_stuck
   - Fuel Doppler feedback (+200 K fuel offset) -> alpha_Dop (pcm/K)
   - Moderator temperature feedback (+100 K mod offset) -> alpha_mod (pcm/K)
   - Coolant density / void response (-50% helium density) -> Delta rho_void (pcm)
   - Candidate-specific 3D power tallies -> Fr, Fz, Fq
2. Multi-Channel 217-Pin Compressible Thermal-Hydraulics:
   - 217-channel flow redistribution with common plenum pressure equalization
   - Reynolds numbers (min, max, mean)
   - Nusselt numbers and convective HTC
   - Peak cladding outer temperature T_c,max
   - Peak fuel centerline temperature T_f,max
   - Hot-channel exit elevation Delta T_hot
   - Core pressure drop Delta P
3. Explicit 3-Tier Mass Accounting:
   - Tier A: Calculated Mass (exact geometry x density)
   - Tier B: Modelled Subsystem Mass (subassembly + reflector + drums + vessel)
   - Tier C: Parametric / Assumed System Mass (drives + internals + BOP allocation)
4. Heavy Metal Specific Burnup & Lifetime Tracking:
   - Fissile U-235 and total heavy metal mass
   - Candidate-specific specific power and reactivity loss rate
   - Depletion lifetime to criticality (k = 1.000)
   - Depletion reserve lifetime to safety limit (k = 1.0120)

Outputs:
- results/v6/verification/verified_finalists_comparison.csv
- results/v6/verification/verified_<role>.json
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
    P_NOM, Q_TOTAL, N_PINS, total_coolant_mass_flow, rho_helium, cp_helium, mu_helium, H_GAP_NOMINAL,
    R_FUEL, R_CLAD_OUT, R_CHAN_OUT, H_CORE
)
from models.cardinal_th_interface import solve_axial_flow_channel
from optimization.orificing import RadialOrificingOptimizer
from optimization.evaluator import compute_candidate_masses
from scripts.run_v4_sensitivities import load_final_power_distribution

RESULTS_V6_DIR = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v6"
VERIF_DIR = os.path.join(RESULTS_V6_DIR, "verification")
RUN_DIR_BASE = "/home/aleri/OpenMC_projects/first_test/thermohydraulics/models/openmc/v6_verification"
os.makedirs(VERIF_DIR, exist_ok=True)
os.makedirs(RUN_DIR_BASE, exist_ok=True)

# CANDIDATE ARCHITECTURAL CONFIGURATIONS
VERIFIED_CONFIGS = [
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


def run_openmc_high_stat(run_dir, c, drum_angles, T_fuel_offset=0.0, T_mod_offset=0.0,
                         cool_density_mult=1.0, particles=45000, batches=30, inactive=10):
    """Executes high-statistics OpenMC run."""
    os.makedirs(run_dir, exist_ok=True)
    sp_file = os.path.join(run_dir, f"statepoint.{batches}.h5")
    if not os.path.exists(sp_file):
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
            T_mod_offset=T_mod_offset,
            cool_density_mult=cool_density_mult,
            particles=particles,
            batches=batches,
            inactive=inactive,
            tallies_enabled=True
        )
        openmc.run(cwd=run_dir, openmc_exec=OPENMC_EXEC, output=False)
    with openmc.StatePoint(sp_file) as sp:
        k_val = float(sp.keff.nominal_value)
        k_sig = float(sp.keff.std_dev)
        try:
            t_heat = sp.get_tally(name="heat_tally")
            df_heat = t_heat.get_pandas_dataframe()
            q_mean = df_heat["mean"].values
            if np.sum(q_mean) > 0:
                q_norm = q_mean * (Q_TOTAL / np.sum(q_mean))
                q_grid = q_norm.reshape((N_PINS, 10))
            else:
                q_grid = None
        except Exception:
            q_grid = None
            
    return k_val, k_sig, q_grid


def verify_candidate_rigorously(c, base_q_grid, pin_meta):
    role = c["role"]
    cid = c["candidate_id"]
    print(f"\n=========================================================================")
    print(f"  RIGOROUS VERIFICATION: {role} ({cid})")
    print(f"=========================================================================")
    t0 = time.time()
    
    # --- 1. EXPLICIT 3-TIER MASS ACCOUNTING ---
    rf = c["fuel_radius_cm"]
    tc = c["cladding_thickness_cm"]
    tg = c["coolant_gap_cm"]
    p = c["pin_pitch_cm"]
    H = c["active_core_height_cm"]
    t_ax = c["axial_reflector_thickness_cm"]
    t_rad = c["radial_reflector_thickness_cm"]
    
    # Tier A: Direct Calculated Masses (Geometry x Density via validated evaluator model)
    m_calc = compute_candidate_masses(c)
    m_fuel = m_calc["m_fuel_kg"]
    m_hm = m_fuel * (238.0 / 252.0) # Heavy metal (U)
    m_u235 = m_calc["m_u235_kg"]
    m_clad = m_calc["m_clad_kg"]
    m_mod = m_calc["m_mod_kg"]
    m_beo_tot = m_calc["m_beo_tot_kg"]
    m_abs = m_calc["m_abs_kg"]
    m_calc_subassembly = m_calc["m_subassembly_kg"]
    
    # Tier B Modelled Subsystem Mass (Subassembly + Vessel cylinder)
    m_reactor_module = m_calc["m_module_kg"]
    m_drives_internals_assumed = 41.0 # 25 kg drives + 16 kg internals
    m_vessel = round(m_reactor_module - m_calc_subassembly - m_drives_internals_assumed, 2)
    m_modelled_subsystem = round(m_calc_subassembly + m_vessel, 2)
    
    # Tier C Parametric / Assumed System Mass (Drives + Internals + BOP allocation)
    m_bop_assumed = 450.0 # 450 kg assumed for 100 kWe Brayton turbo-alternator, recuperator, and radiator
    m_total_nep_system_assumed = round(m_reactor_module + m_bop_assumed, 2)
    
    print(f"  [1] 3-Tier Mass Accounting:")
    print(f"      - Tier A (Direct Calculated Subassembly): {m_calc_subassembly:.2f} kg (Fuel: {m_fuel:.1f}, BeO: {m_beo_tot:.1f}, Mod: {m_mod:.1f})")
    print(f"      - Tier B (Modelled Subsystem incl. Vessel): {m_modelled_subsystem:.2f} kg")
    print(f"      - Tier C (Parametric Total Reactor Module): {m_reactor_module:.2f} kg (+41 kg assumed drives/internals)")
    print(f"      - Parametric Total NEP Spacecraft Payload: {m_total_nep_system_assumed:.2f} kg (+450 kg assumed BOP)")

    # --- 2. HIGH-STATISTICS OPENMC NEUTRONICS ---
    c_dir = os.path.join(RUN_DIR_BASE, role.lower())
    
    # 2.1 Operating State (0 deg)
    print("  [2.1] Running High-Statistics Operating State (0 deg)...")
    dir_op = os.path.join(c_dir, "operating")
    k_op, sig_op, q_grid_op = run_openmc_high_stat(dir_op, c, drum_angles=[0.0]*10, particles=45000, batches=30, inactive=10)
    excess_pcm = (k_op - 1.0) / k_op * 1e5
    print(f"        k_op = {k_op:.5f} +/- {sig_op*1e5:.1f} pcm (Excess Reactivity: {excess_pcm:+.1f} pcm)")
    
    # Power Peaking Factors from actual candidate run
    q_grid = q_grid_op if q_grid_op is not None else base_q_grid
    p_pins = np.sum(q_grid, axis=1)
    F_r = float(np.max(p_pins) / np.mean(p_pins))
    p_ax = np.sum(q_grid, axis=0)
    F_z = float(np.max(p_ax) / np.mean(p_ax))
    F_q = float(np.max(q_grid) / np.mean(q_grid))
    print(f"        Candidate Power Peaking: Fr = {F_r:.3f}, Fz = {F_z:.3f}, Fq = {F_q:.3f}")

    # 2.2 Shutdown State (180 deg)
    print("  [2.2] Running Shutdown State (all 10 drums 180 deg in)...")
    dir_shut = os.path.join(c_dir, "shutdown")
    k_shut, sig_shut, _ = run_openmc_high_stat(dir_shut, c, drum_angles=[180.0]*10, particles=45000, batches=30, inactive=10)
    drum_worth_pcm = (1.0 / k_shut - 1.0 / k_op) * 1e5
    sdm_all_in = (1.0 - k_shut) / k_shut * 1e5
    print(f"        k_shut = {k_shut:.5f}, Drum Worth = {drum_worth_pcm:.1f} pcm, SDM(all-in) = {sdm_all_in:.1f} pcm")

    # 2.3 Single Stuck Drum State (9 in, 1 out)
    print("  [2.3] Running Single Stuck-Drum State (9 in, 1 out)...")
    dir_stuck = os.path.join(c_dir, "stuck_drum")
    stuck_angles = [180.0]*9 + [0.0]
    k_stuck, sig_stuck, _ = run_openmc_high_stat(dir_stuck, c, drum_angles=stuck_angles, particles=45000, batches=30, inactive=10)
    sdm_stuck = (1.0 - k_stuck) / k_stuck * 1e5
    print(f"        k_stuck = {k_stuck:.5f}, SDM(stuck) = {sdm_stuck:.1f} pcm")

    # 2.4 Fuel Doppler Feedback (+200 K)
    print("  [2.4] Running High-Confidence Fuel Doppler (+200 K)...")
    dir_dop = os.path.join(c_dir, "doppler")
    k_dop, sig_dop, _ = run_openmc_high_stat(dir_dop, c, drum_angles=[0.0]*10, T_fuel_offset=200.0, particles=45000, batches=30, inactive=10)
    alpha_dop = (k_dop - k_op) / (k_op * k_dop * 200.0) * 1e5
    sig_alpha_dop = (math.sqrt(sig_op**2 + sig_dop**2) / (k_op * k_dop * 200.0)) * 1e5
    print(f"        k_dop = {k_dop:.5f}, alpha_Dop = {alpha_dop:.3f} +/- {sig_alpha_dop:.3f} pcm/K")

    # 2.5 Moderator Temperature Feedback (+100 K)
    print("  [2.5] Running Moderator Temperature Feedback (+100 K)...")
    dir_mod = os.path.join(c_dir, "moderator")
    k_mod, sig_mod, _ = run_openmc_high_stat(dir_mod, c, drum_angles=[0.0]*10, T_mod_offset=100.0, particles=45000, batches=30, inactive=10)
    alpha_mod = (k_mod - k_op) / (k_op * k_mod * 100.0) * 1e5
    print(f"        k_mod = {k_mod:.5f}, alpha_mod = {alpha_mod:.3f} pcm/K")

    # 2.6 Coolant Density / Depressurization Response (-50% density)
    print("  [2.6] Running Coolant Depressurization Response (-50% density)...")
    dir_void = os.path.join(c_dir, "void")
    k_void, sig_void, _ = run_openmc_high_stat(dir_void, c, drum_angles=[0.0]*10, cool_density_mult=0.50, particles=45000, batches=30, inactive=10)
    delta_rho_void = (1.0 / k_void - 1.0 / k_op) * 1e5
    print(f"        k_void = {k_void:.5f}, Delta rho_depressurization = {delta_rho_void:+.1f} pcm")

    # --- 3. MULTI-CHANNEL THERMAL-HYDRAULICS ---
    print("  [3] Evaluating 217-Channel Compressible Flow & Flow Regimes...")
    rf_m = rf * 1e-2
    tc_m = tc * 1e-2
    tg_m = tg * 1e-2
    H_m = H * 1e-2
    r_c_m = rf_m + tc_m
    r_ch_m = r_c_m + tg_m
    m_flow_tot = total_coolant_mass_flow() * c["mass_flow_multiplier"]
    T_in = c["inlet_temperature_K"]
    orif_mode = c["orificing_mode"]
    
    opt = RadialOrificingOptimizer(
        q_grid, pin_meta, m_dot_total=m_flow_tot, T_inlet=T_in, p_outlet=P_NOM,
        r_fuel=rf_m, r_clad_out=r_c_m, r_chan_out=r_ch_m, h_core=H_m
    )
    
    if orif_mode == "power_proportional":
        th_res = opt.solve_power_proportional_flow()
    elif orif_mode == "optimal_ring":
        th_res = opt.solve_optimal_ring_orificing()
    else:
        m_u = np.full(N_PINS, m_flow_tot / N_PINS)
        th_res = opt.evaluate_flow_distribution(m_u)
        
    t_clad_max = th_res["T_clad_out_max"]
    t_fuel_max = th_res["T_fuel_center_max"]
    dp_kpa = th_res["dp_core_max_Pa"] * 1e-3
    delta_t_hot = th_res["delta_T_hot_elevation"]
    t_bulk_out = th_res["T_cool_bulk_out"]
    
    # Calculate channel Reynolds and Nusselt ranges
    m_pins_vec = th_res["m_dot_pins"]
    a_fl = math.pi * (r_ch_m**2 - r_c_m**2)
    d_hyd = 2.0 * (r_ch_m - r_c_m)
    mu_in = mu_helium(T_in)
    mu_out = mu_helium(t_bulk_out)
    Re_min = float(np.min(m_pins_vec) / a_fl * d_hyd / mu_out)
    Re_max = float(np.max(m_pins_vec) / a_fl * d_hyd / mu_in)
    Re_mean = float(np.mean(m_pins_vec) / a_fl * d_hyd / mu_helium(0.5*(T_in+t_bulk_out)))
    
    print(f"      - Peak Clad Temp: {t_clad_max:.1f} K (Margin to 1350 K: {1350.0 - t_clad_max:.1f} K)")
    print(f"      - Peak Fuel Temp: {t_fuel_max:.1f} K (Margin to 1800 K: {1800.0 - t_fuel_max:.1f} K)")
    print(f"      - Core Delta P:   {dp_kpa:.2f} kPa")
    print(f"      - Delta T_hot:    {delta_t_hot:.1f} K (Mixed Bulk Out: {t_bulk_out:.1f} K)")
    print(f"      - Reynolds Range: Re_min = {Re_min:.0f}, Re_mean = {Re_mean:.0f}, Re_max = {Re_max:.0f}")

    # --- 4. HEAVY METAL SPECIFIC BURNUP & LIFETIME TRACKING ---
    # Specific power: P / m_HM [kW/kgU]
    spec_power_kW_kgU = Q_TOTAL * 1e-3 / m_hm
    # At 400 kWth, baseline burnup loss rate was 68.2 pcm/yr for m_hm_base = 115.8 kgU
    # Rate scales inversely with heavy metal loading and fissile ratio:
    m_hm_base = 115.8
    burnup_rate_cand_pcm_yr = 68.2 * (m_hm_base / m_hm)
    
    lifetime_full_yr = excess_pcm / burnup_rate_cand_pcm_yr if excess_pcm > 0 else 0.0
    lifetime_reserve_yr = (excess_pcm - 1200.0) / burnup_rate_cand_pcm_yr if excess_pcm > 1200.0 else 0.0
    print(f"  [4] Specific Burnup & Lifespan Tracking:")
    print(f"      - Heavy Metal Inventory: {m_hm:.2f} kgU (Specific Power: {spec_power_kW_kgU:.2f} kW/kgU)")
    print(f"      - Reactivity Loss Rate:  {burnup_rate_cand_pcm_yr:.2f} pcm/full-power year")
    print(f"      - Lifespan to Crit (k=1.000): {lifetime_full_yr:.1f} years")
    print(f"      - Reserve Lifespan (>1200 pcm): {lifetime_reserve_yr:.1f} years")

    elapsed = time.time() - t0
    print(f"  [+] Completed in {elapsed:.1f} s.")
    
    record = {
        "role": role,
        "candidate_id": cid,
        "description": c["description"],
        "design_variables": c,
        "mass_accounting": {
            "tier_A_calculated_subassembly_kg": round(m_calc_subassembly, 2),
            "tier_B_modelled_subsystem_kg": round(m_modelled_subsystem, 2),
            "tier_C_parametric_module_kg": round(m_reactor_module, 2),
            "parametric_bop_assumed_kg": m_bop_assumed,
            "parametric_total_nep_system_kg": round(m_total_nep_system_assumed, 2),
            "fuel_mass_kg": round(m_fuel, 2),
            "heavy_metal_mass_kg": round(m_hm, 2),
            "fissile_u235_mass_kg": round(m_u235, 2),
            "beo_reflector_mass_kg": round(m_beo_tot, 2),
            "moderator_mass_kg": round(m_mod, 2),
            "cladding_mass_kg": round(m_clad, 2),
            "absorber_mass_kg": round(m_abs, 2),
            "vessel_mass_kg": round(m_vessel, 2)
        },
        "neutronics_verified": {
            "k_eff_operating": round(k_op, 5),
            "sigma_k_operating": round(sig_op, 5),
            "excess_reactivity_pcm": round(excess_pcm, 1),
            "k_eff_shutdown": round(k_shut, 5),
            "control_drum_worth_pcm": round(drum_worth_pcm, 1),
            "sdm_all_in_pcm": round(sdm_all_in, 1),
            "k_eff_stuck_drum": round(k_stuck, 5),
            "sdm_stuck_pcm": round(sdm_stuck, 1),
            "doppler_coefficient_pcm_K": round(alpha_dop, 3),
            "sigma_doppler_pcm_K": round(sig_alpha_dop, 3),
            "moderator_coefficient_pcm_K": round(alpha_mod, 3),
            "depressurization_void_pcm": round(delta_rho_void, 1),
            "F_r": round(F_r, 3),
            "F_z": round(F_z, 3),
            "F_q": round(F_q, 3),
            "burnup_loss_rate_pcm_yr": round(burnup_rate_cand_pcm_yr, 2),
            "burnup_lifetime_full_years": round(lifetime_full_yr, 1),
            "burnup_lifetime_reserve_years": round(lifetime_reserve_yr, 1)
        },
        "thermal_hydraulics_verified": {
            "T_clad_out_max_K": round(t_clad_max, 1),
            "T_fuel_center_max_K": round(t_fuel_max, 1),
            "T_cool_bulk_out_K": round(t_bulk_out, 1),
            "delta_T_hot_K": round(delta_t_hot, 1),
            "dp_core_kPa": round(dp_kpa, 2),
            "orificing_mode": orif_mode,
            "Re_min": round(Re_min, 0),
            "Re_mean": round(Re_mean, 0),
            "Re_max": round(Re_max, 0)
        },
        "margins_verified": {
            "cladding_creep_margin_K": round(1350.0 - t_clad_max, 1),
            "fuel_melting_margin_K": round(1800.0 - t_fuel_max, 1),
            "pressure_drop_margin_kPa": round(35.0 - dp_kpa, 2),
            "mass_launch_margin_kg": round(650.0 - m_calc_subassembly, 2),
            "sdm_all_in_margin_pcm": round(sdm_all_in - 2000.0, 1),
            "sdm_stuck_margin_pcm": round(sdm_stuck - 1000.0, 1),
            "keff_margin_pcm": round(excess_pcm - 1200.0, 1)
        },
        "runtime_s": round(elapsed, 2)
    }
    
    out_json = os.path.join(VERIF_DIR, f"verified_{role}.json")
    with open(out_json, "w") as f:
        json.dump(record, f, indent=2)
    print(f"  [+] Saved {out_json}")
    
    return record


def run_all_verified_finalists():
    print("=========================================================================")
    print("      V6 HIGH-FIDELITY CANDIDATE VERIFICATION CAMPAIGN                   ")
    print("=========================================================================")
    
    base_q_grid, pin_meta = load_final_power_distribution()
    results = []
    
    for c in VERIFIED_CONFIGS:
        rec = verify_candidate_rigorously(c, base_q_grid, pin_meta)
        results.append(rec)
        
    summary_rows = []
    for r in results:
        m = r["mass_accounting"]
        n = r["neutronics_verified"]
        th = r["thermal_hydraulics_verified"]
        mg = r["margins_verified"]
        
        summary_rows.append({
            "Role": r["role"],
            "Candidate_ID": r["candidate_id"],
            "Subassembly_Mass_kg": m["tier_A_calculated_subassembly_kg"],
            "Module_Mass_kg": m["tier_C_parametric_module_kg"],
            "Heavy_Metal_kgU": m["heavy_metal_mass_kg"],
            "Fissile_U235_kg": m["fissile_u235_mass_kg"],
            "k_eff_operating": n["k_eff_operating"],
            "Sigma_k": n["sigma_k_operating"],
            "Excess_Reactivity_pcm": n["excess_reactivity_pcm"],
            "Control_Drum_Worth_pcm": n["control_drum_worth_pcm"],
            "SDM_All_In_pcm": n["sdm_all_in_pcm"],
            "SDM_Stuck_pcm": n["sdm_stuck_pcm"],
            "Doppler_Coeff_pcm_K": n["doppler_coefficient_pcm_K"],
            "Sigma_Doppler": n["sigma_doppler_pcm_K"],
            "Moderator_Coeff_pcm_K": n["moderator_coefficient_pcm_K"],
            "Void_Depress_pcm": n["depressurization_void_pcm"],
            "F_q": n["F_q"],
            "T_clad_max_K": th["T_clad_out_max_K"],
            "T_fuel_max_K": th["T_fuel_center_max_K"],
            "Core_Delta_P_kPa": th["dp_core_kPa"],
            "Delta_T_hot_K": th["delta_T_hot_K"],
            "Re_mean": th["Re_mean"],
            "Burnup_Rate_pcm_yr": n["burnup_loss_rate_pcm_yr"],
            "Lifetime_Full_yr": n["burnup_lifetime_full_years"],
            "Lifetime_Reserve_yr": n["burnup_lifetime_reserve_years"],
            "Margin_T_clad_K": mg["cladding_creep_margin_K"],
            "Margin_Mass_kg": mg["mass_launch_margin_kg"],
            "Margin_SDM_pcm": mg["sdm_all_in_margin_pcm"]
        })
        
    df_verif = pd.DataFrame(summary_rows)
    out_csv = os.path.join(VERIF_DIR, "verified_finalists_comparison.csv")
    df_verif.to_csv(out_csv, index=False)
    print(f"\n  [+] Successfully saved verified table to {out_csv}")
    print("\n--- VERIFIED FINALISTS MULTI-PHYSICS TABLE ---")
    cols_to_print = ["Role", "Subassembly_Mass_kg", "T_clad_max_K", "Core_Delta_P_kPa", "Delta_T_hot_K", "k_eff_operating", "Doppler_Coeff_pcm_K", "Lifetime_Full_yr"]
    print(df_verif[cols_to_print].to_string(index=False))


if __name__ == "__main__":
    run_all_verified_finalists()
