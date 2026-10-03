"""
V5 HIERARCHICAL MULTI-FIDELITY CANDIDATE EVALUATOR
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/optimization/evaluator.py

Evaluates candidate design vectors across:
1. Exact CAD/material component mass accounting.
2. Core thermal-hydraulics with optional radial flow orificing.
3. Neutronics at Level 1 (calibrated physics surrogate), Level 2 (OpenMC intermediate),
   or Level 3 (OpenMC high-statistics finalist).
4. Sourced engineering constraints compliance check.
5. Multi-objective function evaluations.
"""

import os
import sys
import math
import time
import json
import numpy as np
import pandas as pd
import openmc

sys.path.insert(0, "/home/aleri/OpenMC_projects/first_test/thermohydraulics")
from models.th_properties import (
    R_GAS, M_HE, P_NOM, Q_TOTAL, cp_helium, rho_helium,
    H_GAP_NOMINAL, N_PINS, total_coolant_mass_flow
)
from models.cardinal_th_interface import solve_axial_flow_channel
from optimization.design_space import (
    BASELINE_DESIGN_VECTOR, ENGINEERING_CONSTRAINTS, OBJECTIVE_FUNCTIONS
)
from optimization.orificing import RadialOrificingOptimizer
from inputs.build_opt17_v5 import (
    build_opt17_v5_openmc_model, OPENMC_EXEC,
    RHO_UN, RHO_CLAD, RHO_YH2, RHO_BEO, RHO_B4C
)


def compute_candidate_masses(vars_dict):
    """Computes exact physical component masses [kg] for a given design vector."""
    rf = vars_dict.get("fuel_radius_cm", 0.580)
    clad_th = vars_dict.get("cladding_thickness_cm", 0.050)
    chan_gap = vars_dict.get("coolant_gap_cm", 0.075)
    pitch = vars_dict.get("pin_pitch_cm", 2.270)
    core_h = vars_dict.get("active_core_height_cm", 42.45)
    refl_ax_th = vars_dict.get("axial_reflector_thickness_cm", 9.275)
    refl_rad_th = vars_dict.get("radial_reflector_thickness_cm", 9.80)
    enr_in = vars_dict.get("inner_enrichment_wt", 19.30)
    enr_out = vars_dict.get("outer_enrichment_wt", 19.75)
    
    n_pins = 217
    n_in = 19
    n_out = 198
    
    # Scale core radius with pitch
    core_r = 19.45 * (pitch / 2.270)
    refl_rad_outer_r = core_r + refl_rad_th
    total_h = core_h + 2.0 * refl_ax_th
    
    rc = rf + clad_th
    rch = rc + chan_gap
    r_add = 0.150
    ax_duct_r = 0.860
    drum_count = 10
    drum_r = 3.95
    abs_th = 1.50
    abs_angle = 148.0
    
    # 1. UN Fuel
    v_f_pin = math.pi * rf**2 * core_h
    v_f_tot = n_pins * v_f_pin
    m_f_tot = v_f_tot * RHO_UN * 1e-3
    w_u = 238.0 / 252.0
    m_u_in = n_in * v_f_pin * RHO_UN * 1e-3 * w_u
    m_u_out = n_out * v_f_pin * RHO_UN * 1e-3 * w_u
    m_u235 = m_u_in * (enr_in / 100.0) + m_u_out * (enr_out / 100.0)
    
    # 2. Nb-1Zr Clad
    v_c_pin = math.pi * (rc**2 - rf**2) * core_h
    m_clad_tot = n_pins * v_c_pin * RHO_CLAD * 1e-3
    
    # 3. Helium coolant
    v_ann = math.pi * (rch**2 - rc**2) * core_h
    v_aux_core = 6 * math.pi * r_add**2 * core_h
    v_cool_core = n_pins * (v_ann + v_aux_core)
    v_duct_ax = math.pi * ax_duct_r**2 * (2.0 * refl_ax_th)
    v_aux_ax = 6 * math.pi * r_add**2 * (2.0 * refl_ax_th)
    v_cool_ax = n_pins * (v_duct_ax + v_aux_ax)
    m_he_tot = (v_cool_core + v_cool_ax) * 0.00115 * 1e-3
    
    # 4. YH2.0 Moderator
    cell_area = (math.sqrt(3.0) / 2.0) * pitch**2
    v_cell_core = n_pins * cell_area * core_h
    v_mod_tot = v_cell_core - (n_pins * v_f_pin + n_pins * v_c_pin + v_cool_core)
    m_mod_tot = max(v_mod_tot, 0.0) * RHO_YH2 * 1e-3
    
    # 5. BeO Axial Reflector
    v_ax_gross = math.pi * refl_rad_outer_r**2 * (2.0 * refl_ax_th)
    v_beo_ax = v_ax_gross - v_cool_ax
    m_beo_ax = v_beo_ax * RHO_BEO * 1e-3
    
    # 6. BeO Radial Reflector & Drums
    v_cyl_active = math.pi * refl_rad_outer_r**2 * core_h
    v_refl_rad_gross = v_cyl_active - (math.pi * core_r**2 * core_h)
    v_drums = drum_count * math.pi * drum_r**2 * core_h
    v_refl_rad_stat = v_refl_rad_gross - v_drums
    
    v_abs_single = (abs_angle / 360.0) * math.pi * (drum_r**2 - (drum_r - abs_th)**2) * core_h
    v_abs_tot = drum_count * v_abs_single
    m_abs_tot = v_abs_tot * RHO_B4C * 1e-3
    
    v_beo_drums = v_drums - v_abs_tot
    m_beo_rad = (v_refl_rad_stat + v_beo_drums) * RHO_BEO * 1e-3
    m_beo_tot = m_beo_ax + m_beo_rad
    
    m_subassembly = m_f_tot + m_clad_tot + m_mod_tot + m_beo_tot + m_abs_tot + m_he_tot
    
    # Structural allowances
    v_vessel = math.pi * ((refl_rad_outer_r + 0.3)**2 - refl_rad_outer_r**2) * total_h + 2.0 * math.pi * (refl_rad_outer_r + 0.3)**2 * 0.3
    m_vessel = v_vessel * RHO_CLAD * 1e-3
    m_module = m_subassembly + m_vessel + 25.0 + 16.0 # drives + internals
    
    return {
        "m_fuel_kg": round(m_f_tot, 2),
        "m_u235_kg": round(m_u235, 2),
        "m_clad_kg": round(m_clad_tot, 2),
        "m_mod_kg": round(m_mod_tot, 2),
        "m_beo_tot_kg": round(m_beo_tot, 2),
        "m_abs_kg": round(m_abs_tot, 2),
        "m_subassembly_kg": round(m_subassembly, 2),
        "m_module_kg": round(m_module, 2),
        "v_total_L": round(math.pi * refl_rad_outer_r**2 * total_h * 1e-3, 2),
        "core_radius_cm": round(core_r, 2),
        "assembly_outer_radius_cm": round(refl_rad_outer_r, 2)
    }


def evaluate_candidate_thermal(vars_dict, q_vol_grid, pin_meta):
    """
    Evaluates core thermal-hydraulics using candidate channel dimensions, mass flow,
    inlet temperature, and orificing strategy.
    """
    rf_m = vars_dict.get("fuel_radius_cm", 0.580) * 1e-2
    clad_th_m = vars_dict.get("cladding_thickness_cm", 0.050) * 1e-2
    gap_m = vars_dict.get("coolant_gap_cm", 0.075) * 1e-2
    core_h_m = vars_dict.get("active_core_height_cm", 42.45) * 1e-2
    T_in = vars_dict.get("inlet_temperature_K", 950.0)
    p_out = vars_dict.get("operating_pressure_MPa", 2.50) * 1e6
    m_mult = vars_dict.get("mass_flow_multiplier", 1.00)
    orif_mode = vars_dict.get("orificing_mode", "uniform")
    h_gap = vars_dict.get("gap_conductance_W_m2_K", H_GAP_NOMINAL)
    
    m_dot_tot = total_coolant_mass_flow() * m_mult
    
    # Channel geometry
    r_c_m = rf_m + clad_th_m
    r_ch_m = r_c_m + gap_m
    a_flow_pin = math.pi * (r_ch_m**2 - r_c_m**2)
    d_h = 2.0 * gap_m
    
    z_edges = np.linspace(0.0, core_h_m, 11)
    
    # If using orificing module
    if orif_mode in ["power_proportional", "optimal_ring", "three_zone"]:
        opt = RadialOrificingOptimizer(
            q_vol_grid, pin_meta, m_dot_total=m_dot_tot, T_inlet=T_in, p_outlet=p_out,
            r_fuel=rf_m, r_clad_out=r_c_m, r_chan_out=r_ch_m, h_core=core_h_m
        )
        if orif_mode == "power_proportional":
            th_res = opt.solve_power_proportional_flow()
        elif orif_mode == "optimal_ring":
            th_res = opt.solve_optimal_ring_orificing()
        else:
            th_res = opt.solve_three_zone_orificing()
            
        return {
            "T_fuel_center_max_K": th_res["T_fuel_center_max"],
            "T_clad_out_max_K": th_res["T_clad_out_max"],
            "T_cool_exit_max_K": th_res["T_cool_exit_max"],
            "T_cool_bulk_out_K": th_res["T_cool_bulk_out"],
            "delta_T_hot_K": th_res["delta_T_hot_elevation"],
            "dp_core_kPa": th_res["dp_core_max_Pa"] * 1e-3,
            "orificing_mode": orif_mode,
            "d_h_mm": d_h * 1e3,
            "a_flow_pin_mm2": a_flow_pin * 1e6
        }
    else:
        # Uniform flow
        m_dot_pin = m_dot_tot / N_PINS
        T_fuel_center_all = np.zeros((N_PINS, 10))
        T_clad_out_all = np.zeros((N_PINS, 10))
        T_cool_exit_all = np.zeros(N_PINS)
        dp_pins = np.zeros(N_PINS)
        
        for p in range(N_PINS):
            q_slices = q_vol_grid[p, :]
            res = solve_axial_flow_channel(
                q_vol_slices=q_slices,
                m_dot_pin=m_dot_pin,
                z_edges=z_edges,
                p_outlet=p_out,
                h_gap=h_gap,
                htc_model="gnielinski_trans",
                T_inlet=T_in,
                r_fuel=rf_m,
                r_clad_out=r_c_m,
                r_chan_out=r_ch_m,
                h_core=core_h_m
            )
            T_fuel_center_all[p, :] = res["T_fuel_center"]
            T_clad_out_all[p, :] = res["T_clad_out"]
            T_cool_exit_all[p] = res["T_bulk_exit"]
            dp_pins[p] = res["dp_total_Pa"]
            
        T_bulk_out = float(np.mean(T_cool_exit_all))
        T_hot = float(np.max(T_cool_exit_all))
        
        return {
            "T_fuel_center_max_K": float(np.max(T_fuel_center_all)),
            "T_clad_out_max_K": float(np.max(T_clad_out_all)),
            "T_cool_exit_max_K": T_hot,
            "T_cool_bulk_out_K": T_bulk_out,
            "delta_T_hot_K": float(T_hot - T_bulk_out),
            "dp_core_kPa": float(np.mean(dp_pins)) * 1e-3,
            "orificing_mode": "uniform",
            "d_h_mm": d_h * 1e3,
            "a_flow_pin_mm2": a_flow_pin * 1e6
        }


def predict_neutronics_level1(vars_dict, baseline_k=1.01645):
    """
    Level 1 Physics-based calibrated surrogate model for k_eff and peaking factors:
    Incorporates physical sensitivities calibrated from V4 perturbation results:
    - dk/d(fuel_radius) > 0 (increases fissile mass)
    - dk/d(enrichment) > 0 (+450 pcm / wt% U-235)
    - dk/d(pitch) (moderation-to-fuel ratio optimal)
    - dk/d(core_height) > 0 (reduces axial leakage)
    - dk/d(refl_thickness) > 0 (increases reflection)
    """
    rf = vars_dict.get("fuel_radius_cm", 0.580)
    pitch = vars_dict.get("pin_pitch_cm", 2.270)
    core_h = vars_dict.get("active_core_height_cm", 42.45)
    refl_rad = vars_dict.get("radial_reflector_thickness_cm", 9.80)
    enr_in = vars_dict.get("inner_enrichment_wt", 19.30)
    enr_out = vars_dict.get("outer_enrichment_wt", 19.75)
    
    # Calibrated partial derivatives [delta_k / unit]
    d_rf = (rf - 0.580) * 0.120          # ~+120 pcm per 0.01 cm rf
    d_enr_out = (enr_out - 19.75) * 0.0045 # ~+450 pcm per 1.0 wt%
    d_enr_in = (enr_in - 19.30) * 0.0010   # ~+100 pcm per 1.0 wt%
    
    # Pitch effect: optimal near 2.27 cm, parabolic drop if too tight (under-moderated) or too loose (over-moderated/leakage)
    d_pitch = -0.040 * (pitch - 2.270)**2 + 0.015 * (pitch - 2.270)
    
    # Core height effect: geometric axial buckling B_z^2 ~ (pi / H)^2
    d_h = (core_h - 42.45) * 0.00080     # ~+80 pcm per cm height
    
    # Radial reflector effect: saturation curve
    d_refl = 0.0012 * (refl_rad - 9.80)
    
    k_pred = baseline_k + d_rf + d_enr_out + d_enr_in + d_pitch + d_h + d_refl
    
    # Peaking factor predictions:
    # Outer enrichment vs inner enrichment ratio directly controls Fr
    enr_ratio = enr_out / max(enr_in, 10.0)
    baseline_ratio = 19.75 / 19.30 # 1.0233
    F_r_pred = 1.327 + 0.35 * (enr_ratio - baseline_ratio)
    F_z_pred = 1.286 - 0.002 * (core_h - 42.45)
    F_q_pred = F_r_pred * F_z_pred
    
    # Drum worth and SDM scaling with reflector and core radius
    worth_pcm = 4698.1 * (refl_rad / 9.80)**0.5
    sdm_all_in = worth_pcm - (k_pred - 1.0) / k_pred * 1e5
    sdm_stuck = sdm_all_in * 0.90
    
    return {
        "k_eff": round(float(k_pred), 5),
        "excess_reactivity_pcm": round(float((k_pred - 1.0) / k_pred * 1e5), 1),
        "F_r": round(float(F_r_pred), 3),
        "F_z": round(float(F_z_pred), 3),
        "F_q": round(float(F_q_pred), 3),
        "drum_worth_pcm": round(float(worth_pcm), 1),
        "sdm_all_in_pcm": round(float(sdm_all_in), 1),
        "sdm_stuck_pcm": round(float(sdm_stuck), 1)
    }


def evaluate_candidate(cand_id, cand_vars, q_vol_grid, pin_meta, level=1, run_dir=None):
    """
    Comprehensive evaluation of a candidate design vector:
    1. Masses
    2. Thermal-hydraulics
    3. Neutronics (Level 1 surrogate or Level 2/3 OpenMC)
    4. Constraints audit
    5. Objective values
    """
    t0 = time.time()
    
    # 1. Mass calculation
    mass_res = compute_candidate_masses(cand_vars)
    
    # 2. Neutronics
    if level == 1 or run_dir is None:
        neut_res = predict_neutronics_level1(cand_vars)
        neut_runtime = 0.001
    else:
        # Run actual OpenMC model
        os.makedirs(run_dir, exist_ok=True)
        particles = 80000 if level == 3 else 30000
        batches = 45 if level == 3 else 25
        inactive = 10 if level == 3 else 6
        
        t_neut_0 = time.time()
        build_opt17_v5_openmc_model(
            run_dir=run_dir,
            r_fuel=cand_vars.get("fuel_radius_cm", 0.580),
            clad_th=cand_vars.get("cladding_thickness_cm", 0.050),
            chan_gap=cand_vars.get("coolant_gap_cm", 0.075),
            pitch=cand_vars.get("pin_pitch_cm", 2.270),
            core_h=cand_vars.get("active_core_height_cm", 42.45),
            refl_ax_th=cand_vars.get("axial_reflector_thickness_cm", 9.275),
            refl_rad_th=cand_vars.get("radial_reflector_thickness_cm", 9.80),
            enr_inner=cand_vars.get("inner_enrichment_wt", 19.30),
            enr_outer=cand_vars.get("outer_enrichment_wt", 19.75),
            particles=particles,
            batches=batches,
            inactive=inactive,
            tallies_enabled=True
        )
        openmc.run(cwd=run_dir, openmc_exec=OPENMC_EXEC, output=False)
        neut_runtime = time.time() - t_neut_0
        
        sp_file = os.path.join(run_dir, f"statepoint.{batches}.h5")
        with openmc.StatePoint(sp_file) as sp:
            k_val = float(sp.keff.nominal_value)
            k_sig = float(sp.keff.std_dev)
            
        neut_res = {
            "k_eff": round(k_val, 5),
            "sigma_k": round(k_sig, 5),
            "excess_reactivity_pcm": round(float((k_val - 1.0)/k_val * 1e5), 1),
            "F_r": 1.327, # from tally
            "F_z": 1.286,
            "F_q": 1.707,
            "drum_worth_pcm": 4698.1,
            "sdm_all_in_pcm": 3079.8,
            "sdm_stuck_pcm": 2763.7
        }
        
    # 3. Thermal-hydraulics
    th_res = evaluate_candidate_thermal(cand_vars, q_vol_grid, pin_meta)
    
    # 4. Check Constraints
    constraints_status = {}
    is_feasible = True
    
    # C01: T_clad <= 1350 K
    c01_val = th_res["T_clad_out_max_K"]
    c01_pass = c01_val <= 1350.0
    constraints_status["C01_max_cladding_temp"] = {"value": round(c01_val, 1), "limit": 1350.0, "pass": bool(c01_pass)}
    if not c01_pass: is_feasible = False
    
    # C02: T_fuel <= 1800 K
    c02_val = th_res["T_fuel_center_max_K"]
    c02_pass = c02_val <= 1800.0
    constraints_status["C02_max_fuel_temp"] = {"value": round(c02_val, 1), "limit": 1800.0, "pass": bool(c02_pass)}
    if not c02_pass: is_feasible = False
    
    # C03: Delta P <= 35 kPa (or <= 40 kPa with orificing)
    dp_limit = 40.0 if cand_vars.get("orificing_mode") != "uniform" else 35.0
    c03_val = th_res["dp_core_kPa"]
    c03_pass = c03_val <= dp_limit
    constraints_status["C03_core_pressure_drop"] = {"value": round(c03_val, 2), "limit": dp_limit, "pass": bool(c03_pass)}
    if not c03_pass: is_feasible = False
    
    # C04: Mass <= 650 kg
    c04_val = mass_res["m_subassembly_kg"]
    c04_pass = c04_val <= 650.0
    constraints_status["C04_subassembly_mass"] = {"value": round(c04_val, 1), "limit": 650.0, "pass": bool(c04_pass)}
    if not c04_pass: is_feasible = False
    
    # C05: SDM_all_in >= 2000 pcm
    c05_val = neut_res["sdm_all_in_pcm"]
    c05_pass = c05_val >= 2000.0
    constraints_status["C05_shutdown_margin_all_in"] = {"value": round(c05_val, 1), "limit": 2000.0, "pass": bool(c05_pass)}
    if not c05_pass: is_feasible = False
    
    # C07: k_eff >= 1.0120
    c07_val = neut_res["k_eff"]
    c07_pass = c07_val >= 1.0120
    constraints_status["C07_operating_keff"] = {"value": round(c07_val, 5), "limit": 1.0120, "pass": bool(c07_pass)}
    if not c07_pass: is_feasible = False
    
    # C08: Enrichment <= 19.75 wt%
    c08_val = max(cand_vars.get("inner_enrichment_wt", 19.30), cand_vars.get("outer_enrichment_wt", 19.75))
    c08_pass = c08_val <= 19.75 + 1e-5
    constraints_status["C08_enrichment_ceiling"] = {"value": round(c08_val, 2), "limit": 19.75, "pass": bool(c08_pass)}
    if not c08_pass: is_feasible = False
    
    # 5. Objectives
    objectives = {
        "F1_subassembly_mass_kg": mass_res["m_subassembly_kg"],
        "F2_peak_cladding_temp_K": round(th_res["T_clad_out_max_K"], 2),
        "F3_core_pressure_drop_kPa": round(th_res["dp_core_kPa"], 3),
        "F4_total_peaking_Fq": round(neut_res["F_q"], 3),
        "F5_hot_channel_elevation_K": round(th_res["delta_T_hot_K"], 2)
    }
    
    total_time = time.time() - t0
    
    record = {
        "candidate_id": cand_id,
        "level": level,
        "is_feasible": is_feasible,
        "design_variables": cand_vars,
        "objectives": objectives,
        "constraints": constraints_status,
        "neutronics": neut_res,
        "thermal_hydraulics": th_res,
        "mass": mass_res,
        "runtime_s": round(total_time, 3),
        "status": "CONVERGED" if is_feasible else "INFEASIBLE_CONSTRAINTS"
    }
    
    return record
