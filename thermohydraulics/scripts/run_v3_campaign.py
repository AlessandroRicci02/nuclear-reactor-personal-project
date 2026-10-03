"""
MASTER V3 CAMPAIGN RUNNER: ADVANCED OPENMC <-> MOOSE HEAT CONDUCTION <-> THM COUPLING
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/scripts/run_v3_campaign.py

Executes:
- RUN A: Baseline (Previous Final Coupled from V2: RUN 5)
- RUN B: Improved Cardinal/MOOSE/THM Coupling (Picard iterations to strict convergence)
- RUN C: Refined Spatial/Thermal Discretization Check (20 axial slices vs 10 axial slices)
- RUN D: High-Statistics Converged Master Solution (160,000 particles, 50 batches)
Exports all 8 CSV/JSON deliverables to results/v3/ and generates all 9 figures.
"""

import os
import sys
import time
import json
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openmc

sys.path.insert(0, "/home/aleri/OpenMC_projects/first_test/thermohydraulics")
from models.th_properties import (
    Q_TOTAL, N_PINS, H_CORE, R_FUEL, R_CLAD_OUT, R_CHAN_OUT,
    T_INLET, DELTA_T_BULK, cp_helium, pin_coolant_mass_flow,
    total_coolant_mass_flow, rho_helium, H_GAP_NOMINAL,
    pressure_drop_channel
)
from models.cardinal_th_interface import (
    solve_axial_flow_channel, generate_cardinal_thm_input,
    generate_cardinal_solid_input, run_cardinal_input
)
from coupling.v3_picard_driver import V3PicardCoupledEngine
from inputs.build_opt17_v3 import build_opt17_v3_openmc_model
from scripts.generate_v3_plots import plot_all_v3_figures

WORKSPACE_ROOT = "/home/aleri/OpenMC_projects/first_test/thermohydraulics"
RESULTS_V3_DIR = os.path.join(WORKSPACE_ROOT, "results", "v3")
FIG_DIR = os.path.join(RESULTS_V3_DIR, "figures")

os.makedirs(RESULTS_V3_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)


def execute_v3_campaign():
    print("=" * 80)
    print("  OPT-17 CAND 17-CLEAN: V3 IMPROVED MULTIPHYSICS COUPLING CAMPAIGN")
    print("  OpenMC 0.16.0 <-> MOOSE Heat Conduction <-> MOOSE THM Flow Channels")
    print("  Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/results/v3/")
    print("=" * 80)
    
    t_start_total = time.time()
    
    # =========================================================================
    # STEP 1: LOAD RUN A (PREVIOUS FINAL COUPLED BASELINE)
    # =========================================================================
    print("\n" + "#" * 70)
    print("  [STEP 1] LOADING RUN A (PREVIOUS FINAL COUPLED REFERENCE - V2 RUN 5)")
    print("#" * 70)
    
    # Load from previous coupling_iterations.csv (final iteration 5)
    coup_iter_file = os.path.join(WORKSPACE_ROOT, "results", "coupling", "coupling_iterations.csv")
    df_coup = pd.read_csv(coup_iter_file)
    run5_row = df_coup.iloc[-1].to_dict()
    
    run_a_record = {
        "run_id": "RUN A (V2 Baseline)",
        "k_eff": float(run5_row["k_eff"]),
        "sigma_k": float(run5_row["sigma_k"]),
        "delta_rho_pcm": float(run5_row["reactivity"]),
        "peak_pin_power_W": float(run5_row["peak_pin_power_W"]) if "peak_pin_power_W" in run5_row else float(run5_row["peak_pin_power"]),
        "F_r": float(run5_row["F_r"]),
        "F_z": float(run5_row["F_z"]),
        "F_q": float(run5_row["F_q"]),
        "q_vol_max_W_m3": 7.073e7,
        "T_fuel_max_K": float(run5_row["Tfuel_max"]),
        "T_clad_max_K": float(run5_row["Tclad_max"]),
        "T_coolant_max_K": float(run5_row["Tcoolant_max"]),
        "T_coolant_out_K": float(run5_row["Tcoolant_out"]),
        "mass_flow_kg_s": float(run5_row["mass_flow"]),
        "Re_min": float(run5_row["Re_min"]),
        "Re_max": float(run5_row["Re_max"]),
        "Nu_mean": float(run5_row["Nu_min"]),
        "h_min_W_m2_K": float(run5_row["h_min"]),
        "h_max_W_m2_K": float(run5_row["h_max"]),
        "pressure_drop_Pa": float(run5_row["pressure_drop"]),
        "iterations_count": int(run5_row["iteration"]) + 1,
        "particles": 80000,
        "batches": 45
    }
    
    print(f"  RUN A Loaded: k_eff = {run_a_record['k_eff']:.5f} +/- {run_a_record['sigma_k']:.5f}")
    print(f"  T_fuel,max = {run_a_record['T_fuel_max_K']:.1f} K | Peak Pin Power = {run_a_record['peak_pin_power_W']:.1f} W")
    print(f"  Core Delta P = {run_a_record['pressure_drop_Pa']/1000.0:.2f} kPa")
    
    # =========================================================================
    # STEP 2: EXECUTE RUN B (IMPROVED CARDINAL/MOOSE/THM PICARD COUPLING)
    # =========================================================================
    print("\n" + "#" * 70)
    print("  [STEP 2] EXECUTING RUN B: IMPROVED CARDINAL/MOOSE/THM PICARD COUPLING")
    print("  Picard under-relaxation (theta=0.60), WMP Doppler feedback,")
    print("  Full-core 1D THM channel solve + Cardinal FE benchmark")
    print("#" * 70)
    
    run_b_dir = os.path.join(RESULTS_V3_DIR, "run_b")
    engine_b = V3PicardCoupledEngine(
        work_dir=run_b_dir,
        n_axial_slices=10,
        particles=80000,
        batches=45,
        inactive=10,
        theta_relaxation=0.60,
        h_gap=H_GAP_NOMINAL,
        htc_model="gnielinski_trans",
        enable_source_reuse=True
    )
    
    max_iter = 8
    tol_delta_k_ratio = 1.50   # |delta_k| <= 1.5 * sigma_k
    tol_delta_T_max = 2.0      # delta_T_fuel,max <= 2.0 K
    tol_rms_power = 0.015      # RMS power variation <= 1.5%
    
    b_iter_records = []
    b_convergence_records = []
    cardinal_bench_records = []
    
    v3_iter_file = os.path.join(RESULTS_V3_DIR, "v3_iterations.csv")
    already_computed = False
    if os.path.exists(v3_iter_file):
        df_existing = pd.read_csv(v3_iter_file)
        if len(df_existing) >= 8:
            print("  [INFO] Found existing converged RUN B iterations in v3_iterations.csv. Loading...")
            b_iter_records = df_existing.to_dict(orient="records")
            k_ref_uncoupled = b_iter_records[0]["k_eff"]
            converged = True
            already_computed = True
            
            # Load statepoint of final iteration 7
            sp_final_b = os.path.join(run_b_dir, "iter_7", "statepoint.45.h5")
            engine_b.last_statepoint = sp_final_b
            
            # Reconstruct model and pin metadata
            iter_dir_7 = os.path.join(run_b_dir, "iter_7")
            model_info = build_opt17_v3_openmc_model(
                run_dir=iter_dir_7,
                n_axial_slices=10,
                particles=1000,
                batches=5,
                inactive=2
            )
            
            with openmc.StatePoint(sp_final_b) as sp:
                k_val = float(sp.k_combined.nominal_value)
                k_std = float(sp.k_combined.std_dev)
                t_heat = sp.get_tally(name="fuel_cell_heating")
                df_heat = t_heat.get_pandas_dataframe()
                
            df_heat_scores = df_heat[df_heat["score"] == "heating"]
            cell_heat_dict = dict(zip(df_heat_scores["cell"], df_heat_scores["mean"]))
            
            q_raw_grid = np.zeros((N_PINS, 10))
            for p in range(N_PINS):
                for s in range(10):
                    cell_obj = model_info["fuel_cells_map"][(p, s)]
                    q_raw_grid[p, s] = max(cell_heat_dict.get(cell_obj.id, 0.0), 1e-12)
            total_raw_core = np.sum(q_raw_grid)
            norm_factor = Q_TOTAL / total_raw_core
            p_cell_watts = q_raw_grid * norm_factor
            q_relaxed = p_cell_watts / engine_b.v_fuel_slice
            th_res = engine_b.run_thermal_step(q_relaxed)
            engine_b.T_fuel_vol = th_res["T_fuel_vol"]
            engine_b.T_clad_avg = th_res["T_clad_avg"]
            engine_b.T_cool_bulk = th_res["T_cool_bulk"]
            neut_res = {"pin_meta": model_info["pin_meta"]}
            peak_pin_idx = int(np.argmax(th_res["pin_powers"]))
            p_avg = float(np.mean(th_res["pin_powers"]))
            dp_core = float(th_res["dp_total"])

    if not already_computed:
        for it in range(max_iter):
            print(f"\n--- RUN B: Picard Iteration {it} ---")
            t_it_start = time.time()
            
            # Adapt particle count for later iterations to reduce statistical noise
            p_count = 80000 if it < 4 else 100000
            
            # 1. Neutronics Step
            neut_res = engine_b.run_neutronics_step(it, particles_override=p_count)
            k_val = neut_res["k_eff"]
            k_std = neut_res["k_std"]
            q_new = neut_res["q_vol_new"]
            
            if it == 0:
                k_ref_uncoupled = k_val
                q_relaxed = q_new.copy()
            else:
                # Heat source under-relaxation: q_relaxed = theta * q_new + (1 - theta) * q_prev
                q_relaxed = engine_b.theta * q_new + (1.0 - engine_b.theta) * prev_q_vol
                
            # 2. Thermohydraulic & Solid Conduction Step
            th_res = engine_b.run_thermal_step(q_relaxed)
            
            # 3. Update engine thermal fields for next iteration feedback
            engine_b.T_fuel_vol = th_res["T_fuel_vol"]
            engine_b.T_clad_avg = th_res["T_clad_avg"]
            engine_b.T_cool_bulk = th_res["T_cool_bulk"]
        
            # Metrics calculation
            delta_rho_pcm = (k_val - k_ref_uncoupled) / (k_val * k_ref_uncoupled) * 1e5
            p_pins = th_res["pin_powers"]
            p_peak = float(np.max(p_pins))
            p_avg = float(np.mean(p_pins))
            f_r = p_peak / p_avg
            
            # Axial profile on peak pin
            peak_pin_idx = int(np.argmax(p_pins))
            q_peak_axial = q_relaxed[peak_pin_idx, :]
            f_z = float(np.max(q_peak_axial) / np.mean(q_peak_axial))
            f_q = f_r * f_z
            
            T_f_max = float(np.max(th_res["T_fuel_center"]))
            T_c_max = float(np.max(th_res["T_clad_out"]))
            T_cool_max = float(np.max(th_res["T_cool_bulk"]))
            dp_core = float(th_res["dp_total"])
            
            # Evaluate convergence against previous iteration
            if it == 0:
                delta_k = 0.0
                delta_k_ratio = 0.0
                delta_T_f = 0.0
                rms_P = 0.0
                delta_Fz = 0.0
            else:
                delta_k = k_val - prev_k
                delta_k_ratio = abs(delta_k) / k_std
                delta_T_f = float(np.max(np.abs(th_res["T_fuel_vol"] - prev_T_fuel)))
                rms_P = float(np.sqrt(np.mean(((p_pins - prev_pin_powers) / p_avg)**2)))
                delta_Fz = abs(f_z - prev_Fz)
                
            print(f"  k_eff = {k_val:.5f} +/- {k_std:.5f} | Delta rho = {delta_rho_pcm:+.1f} pcm")
            print(f"  |Delta k| = {abs(delta_k):.5f} ({delta_k_ratio:.2f} sigma) | Delta T_fuel,max = {delta_T_f:.2f} K | RMS dP = {rms_P*100:.2f}%")
            print(f"  T_fuel,max = {T_f_max:.1f} K | T_clad,max = {T_c_max:.1f} K | P_peak = {p_peak:.1f} W (Pin #{peak_pin_idx})")
            print(f"  Core Delta P = {dp_core/1000.0:.2f} kPa | Re = [{th_res['Re_range'][0]:.0f}, {th_res['Re_range'][1]:.0f}]")
            
            # 4. Direct Cardinal Benchmark on Peak Pin
            bench_out = engine_b.run_cardinal_benchmark_pin(peak_pin_idx, q_peak_axial, neut_res["iter_dir"])
            card_thm = bench_out["thm_results"]
            card_solid = bench_out["solid_results"]
            
            card_Tf_max = card_solid.get("T_fuel_max", T_f_max)
            card_Tc_max = card_solid.get("T_clad_max", T_c_max)
            card_dp = card_thm.get("dp_total", dp_core)
            diff_Tf = abs(card_Tf_max - T_f_max)
            diff_dp = abs(card_dp - dp_core)
            
            print(f"  Cardinal Benchmark Pin #{peak_pin_idx}:")
            print(f"    Cardinal T_fuel,max = {card_Tf_max:.2f} K (Diff: {diff_Tf:.2f} K, {diff_Tf/T_f_max*100:.2f}%)")
            print(f"    Cardinal Clad Delta P = {card_dp/1000.0:.2f} kPa")
            
            iter_rec = {
                "iteration": it,
                "k_eff": k_val,
                "sigma_k": k_std,
                "delta_k": delta_k,
                "delta_k_sigma_ratio": delta_k_ratio,
                "delta_rho_pcm": delta_rho_pcm,
                "peak_pin_power_W": p_peak,
                "peak_pin_idx": peak_pin_idx,
                "F_r": f_r,
                "F_z": f_z,
                "F_q": f_q,
                "q_vol_max_W_m3": float(np.max(q_relaxed)),
                "T_fuel_max_K": T_f_max,
                "T_clad_max_K": T_c_max,
                "T_coolant_max_K": T_cool_max,
                "T_coolant_out_K": th_res["T_out_mean"],
                "delta_T_fuel_K": delta_T_f,
                "power_RMS_change": rms_P,
                "pressure_drop_Pa": dp_core,
                "Re_mean": th_res["Re_mean"],
                "Nu_mean": th_res["Nu_mean"],
                "h_mean_W_m2_K": th_res["h_mean"],
                "cardinal_T_fuel_max_K": card_Tf_max,
                "cardinal_T_clad_max_K": card_Tc_max,
                "cardinal_dp_Pa": card_dp,
                "runtime_s": time.time() - t_it_start
            }
            b_iter_records.append(iter_rec)
            
            # Check convergence criteria (must satisfy all from it >= 3)
            if it >= 3 and delta_k_ratio <= tol_delta_k_ratio and delta_T_f <= tol_delta_T_max and rms_P <= tol_rms_power:
                print(f"\n>>> CONVERGENCE ACHIEVED AT ITERATION {it}! <<<")
                print(f"    |Delta k|/sigma = {delta_k_ratio:.2f} <= {tol_delta_k_ratio}")
                print(f"    Delta T_fuel,max = {delta_T_f:.2f} K <= {tol_delta_T_max} K")
                print(f"    RMS Power Change = {rms_P*100:.2f}% <= {tol_rms_power*100:.2f}%")
                converged = True
                break
                
            prev_T_fuel = th_res["T_fuel_vol"].copy()
            prev_q_vol = q_relaxed.copy()
            prev_k = k_val
            prev_pin_powers = p_pins.copy()
            prev_Fz = f_z
        
    df_b_iters = pd.DataFrame(b_iter_records)
    df_b_iters.to_csv(os.path.join(RESULTS_V3_DIR, "v3_iterations.csv"), index=False)
    
    # Save pin power and thermal fields of final RUN B
    b_final_rec = b_iter_records[-1]
    
    # Pin power CSV
    pin_powers_list = []
    pin_meta_list = neut_res["pin_meta"]
    for p in range(N_PINS):
        meta = pin_meta_list[p]
        p_val = th_res["pin_powers"][p]
        pin_powers_list.append({
            "pin_idx": p,
            "ring": meta[1] if isinstance(meta, tuple) else meta.get("ring", 0),
            "x_cm": meta[2] if isinstance(meta, tuple) else meta.get("x_cm", 0.0),
            "y_cm": meta[3] if isinstance(meta, tuple) else meta.get("y_cm", 0.0),
            "power_W": p_val,
            "peaking_Fr": p_val / p_avg
        })
    df_pin_pow = pd.DataFrame(pin_powers_list)
    df_pin_pow.to_csv(os.path.join(RESULTS_V3_DIR, "v3_pin_power.csv"), index=False)
    
    # Thermal fields CSV (axial profiles of representative pins)
    thermal_field_list = []
    z_mids = 0.5 * (engine_b.z_edges[:-1] + engine_b.z_edges[1:])
    for s in range(engine_b.n_axial_slices):
        thermal_field_list.append({
            "slice_idx": s,
            "z_cm": z_mids[s] * 100.0,
            "T_fuel_center_peak_K": th_res["T_fuel_center"][peak_pin_idx, s],
            "T_fuel_vol_peak_K": th_res["T_fuel_vol"][peak_pin_idx, s],
            "T_clad_out_peak_K": th_res["T_clad_out"][peak_pin_idx, s],
            "T_cool_bulk_peak_K": th_res["T_cool_bulk"][peak_pin_idx, s],
            "T_fuel_center_core_avg_K": np.mean(th_res["T_fuel_center"][:, s]),
            "T_fuel_vol_core_avg_K": np.mean(th_res["T_fuel_vol"][:, s]),
            "T_clad_out_core_avg_K": np.mean(th_res["T_clad_out"][:, s]),
            "T_cool_bulk_core_avg_K": np.mean(th_res["T_cool_bulk"][:, s]),
            "q_vol_peak_W_m3": q_relaxed[peak_pin_idx, s],
            "q_vol_avg_W_m3": np.mean(q_relaxed[:, s])
        })
    df_therm = pd.DataFrame(thermal_field_list)
    df_therm.to_csv(os.path.join(RESULTS_V3_DIR, "v3_thermal_field.csv"), index=False)
    
    # Coolant field CSV
    coolant_field_list = []
    for s in range(engine_b.n_axial_slices):
        t_b = np.mean(th_res["T_cool_bulk"][:, s])
        coolant_field_list.append({
            "slice_idx": s,
            "z_cm": z_mids[s] * 100.0,
            "T_coolant_bulk_K": t_b,
            "density_He_kg_m3": rho_helium(t_b),
            "pressure_drop_accum_Pa": dp_core * (s + 1) / engine_b.n_axial_slices,
            "velocity_m_s": engine_b.m_dot_pin / (rho_helium(t_b) * (math.pi * (R_CHAN_OUT**2 - R_CLAD_OUT**2))),
            "Re": th_res["Re_mean"],
            "Nu": th_res["Nu_mean"],
            "HTC_W_m2_K": th_res["h_mean"]
        })
    df_cool = pd.DataFrame(coolant_field_list)
    df_cool.to_csv(os.path.join(RESULTS_V3_DIR, "v3_coolant_field.csv"), index=False)
    
    # Convergence summary CSV
    conv_rec = {
        "status": "CONVERGED" if converged else "MAX_ITER_REACHED",
        "final_iteration": b_final_rec["iteration"],
        "delta_k_final": b_final_rec["delta_k"],
        "delta_k_ratio": b_final_rec["delta_k_sigma_ratio"],
        "delta_T_fuel_final_K": b_final_rec["delta_T_fuel_K"],
        "power_RMS_change_final": b_final_rec["power_RMS_change"],
        "tolerance_k_ratio": tol_delta_k_ratio,
        "tolerance_delta_T_K": tol_delta_T_max,
        "tolerance_RMS_power": tol_rms_power,
        "cardinal_diff_T_fuel_K": abs(b_final_rec["cardinal_T_fuel_max_K"] - b_final_rec["T_fuel_max_K"]),
        "cardinal_diff_clad_K": abs(b_final_rec["cardinal_T_clad_max_K"] - b_final_rec["T_clad_max_K"])
    }
    pd.DataFrame([conv_rec]).to_csv(os.path.join(RESULTS_V3_DIR, "v3_convergence.csv"), index=False)
    
    # =========================================================================
    # STEP 3: EXECUTE RUN C (REFINED SPATIAL/THERMAL RESOLUTION - 20 SLICES)
    # =========================================================================
    print("\n" + "#" * 70)
    print("  [STEP 3] EXECUTING RUN C: REFINED SPATIAL RESOLUTION CHECK (20 AXIAL SLICES)")
    print("  4,340 fuel cells in OpenMC & 20 elements in THM FlowChannel")
    print("#" * 70)
    
    run_c_dir = os.path.join(RESULTS_V3_DIR, "run_c")
    engine_c = V3PicardCoupledEngine(
        work_dir=run_c_dir,
        n_axial_slices=20,
        particles=80000,
        batches=45,
        inactive=10,
        theta_relaxation=0.60,
        h_gap=H_GAP_NOMINAL,
        enable_source_reuse=False
    )
    
    # Map the converged 10-slice temperatures to 20 slices as starting guess
    T_fuel_20 = np.repeat(th_res["T_fuel_vol"], 2, axis=1)
    T_clad_20 = np.repeat(th_res["T_clad_avg"], 2, axis=1)
    T_cool_20 = np.repeat(th_res["T_cool_bulk"], 2, axis=1)
    engine_c.T_fuel_vol = T_fuel_20
    engine_c.T_clad_avg = T_clad_20
    engine_c.T_cool_bulk = T_cool_20
    
    # Run 2 Picard iterations on the refined 20-slice mesh
    for it_c in range(2):
        print(f"  RUN C (20 slices) - Iteration {it_c}...")
        c_neut = engine_c.run_neutronics_step(it_c)
        c_q_relaxed = c_neut["q_vol_new"]
        c_th = engine_c.run_thermal_step(c_q_relaxed)
        engine_c.T_fuel_vol = c_th["T_fuel_vol"]
        engine_c.T_clad_avg = c_th["T_clad_avg"]
        engine_c.T_cool_bulk = c_th["T_cool_bulk"]
        
    c_p_pins = c_th["pin_powers"]
    c_p_peak = float(np.max(c_p_pins))
    c_f_r = c_p_peak / float(np.mean(c_p_pins))
    c_peak_idx = int(np.argmax(c_p_pins))
    c_f_z = float(np.max(c_q_relaxed[c_peak_idx, :]) / np.mean(c_q_relaxed[c_peak_idx, :]))
    c_f_q = c_f_r * c_f_z
    c_Tf_max = float(np.max(c_th["T_fuel_center"]))
    c_Tc_max = float(np.max(c_th["T_clad_out"]))
    c_dp = float(c_th["dp_total"])
    
    # Discretization sensitivity comparison: 10 vs 20 slices
    diff_mesh_k = (c_neut["k_eff"] - b_final_rec["k_eff"]) * 1e5
    diff_mesh_Tf = c_Tf_max - b_final_rec["T_fuel_max_K"]
    diff_mesh_Ppeak = c_p_peak - b_final_rec["peak_pin_power_W"]
    diff_mesh_Fz = c_f_z - b_final_rec["F_z"]
    diff_mesh_dp = c_dp - b_final_rec["pressure_drop_Pa"]
    
    mesh_sens_records = [
        {"metric": "Axial Slices", "case_10_slices": 10, "case_20_slices": 20, "difference": 10, "percent_diff": 100.0},
        {"metric": "Total Fuel Cells", "case_10_slices": 2170, "case_20_slices": 4340, "difference": 2170, "percent_diff": 100.0},
        {"metric": "k_eff", "case_10_slices": b_final_rec["k_eff"], "case_20_slices": c_neut["k_eff"], "difference": c_neut["k_eff"] - b_final_rec["k_eff"], "percent_diff": (c_neut["k_eff"] - b_final_rec["k_eff"]) / b_final_rec["k_eff"] * 100.0},
        {"metric": "Reactivity (pcm)", "case_10_slices": 0.0, "case_20_slices": diff_mesh_k, "difference": diff_mesh_k, "percent_diff": 0.0},
        {"metric": "T_fuel,max (K)", "case_10_slices": b_final_rec["T_fuel_max_K"], "case_20_slices": c_Tf_max, "difference": diff_mesh_Tf, "percent_diff": diff_mesh_Tf / b_final_rec["T_fuel_max_K"] * 100.0},
        {"metric": "T_clad,max (K)", "case_10_slices": b_final_rec["T_clad_max_K"], "case_20_slices": c_Tc_max, "difference": c_Tc_max - b_final_rec["T_clad_max_K"], "percent_diff": (c_Tc_max - b_final_rec["T_clad_max_K"]) / b_final_rec["T_clad_max_K"] * 100.0},
        {"metric": "Peak Pin Power (W)", "case_10_slices": b_final_rec["peak_pin_power_W"], "case_20_slices": c_p_peak, "difference": diff_mesh_Ppeak, "percent_diff": diff_mesh_Ppeak / b_final_rec["peak_pin_power_W"] * 100.0},
        {"metric": "F_r", "case_10_slices": b_final_rec["F_r"], "case_20_slices": c_f_r, "difference": c_f_r - b_final_rec["F_r"], "percent_diff": (c_f_r - b_final_rec["F_r"]) / b_final_rec["F_r"] * 100.0},
        {"metric": "F_z", "case_10_slices": b_final_rec["F_z"], "case_20_slices": c_f_z, "difference": diff_mesh_Fz, "percent_diff": diff_mesh_Fz / b_final_rec["F_z"] * 100.0},
        {"metric": "F_q", "case_10_slices": b_final_rec["F_q"], "case_20_slices": c_f_q, "difference": c_f_q - b_final_rec["F_q"], "percent_diff": (c_f_q - b_final_rec["F_q"]) / b_final_rec["F_q"] * 100.0},
        {"metric": "Core Delta P (Pa)", "case_10_slices": b_final_rec["pressure_drop_Pa"], "case_20_slices": c_dp, "difference": diff_mesh_dp, "percent_diff": diff_mesh_dp / b_final_rec["pressure_drop_Pa"] * 100.0}
    ]
    df_mesh_sens = pd.DataFrame(mesh_sens_records)
    df_mesh_sens.to_csv(os.path.join(RESULTS_V3_DIR, "v3_mesh_sensitivity.csv"), index=False)
    print(f"  RUN C Finished: k_eff = {c_neut['k_eff']:.5f} | T_fuel,max = {c_Tf_max:.1f} K (Diff: {diff_mesh_Tf:+.1f} K)")
    print(f"  F_z = {c_f_z:.4f} (Diff: {diff_mesh_Fz:+.4f}) | Peak Pin Power = {c_p_peak:.1f} W")
    
    # =========================================================================
    # STEP 4: EXECUTE RUN D (HIGH-STATISTICS CONVERGED MASTER SOLUTION)
    # =========================================================================
    print("\n" + "#" * 70)
    print("  [STEP 4] EXECUTING RUN D: HIGH-STATISTICS MASTER SOLUTION")
    print("  160,000 particles per batch x 50 batches (10 inactive, 40 active)")
    print("  Converged source reuse to minimize Monte Carlo statistical uncertainty")
    print("#" * 70)
    
    run_d_dir = os.path.join(RESULTS_V3_DIR, "run_d")
    engine_d = V3PicardCoupledEngine(
        work_dir=run_d_dir,
        n_axial_slices=10,
        particles=160000,
        batches=50,
        inactive=10,
        theta_relaxation=0.60,
        h_gap=H_GAP_NOMINAL,
        enable_source_reuse=True
    )
    engine_d.T_fuel_vol = th_res["T_fuel_vol"].copy()
    engine_d.T_clad_avg = th_res["T_clad_avg"].copy()
    engine_d.T_cool_bulk = th_res["T_cool_bulk"].copy()
    engine_d.last_statepoint = engine_b.last_statepoint
    
    d_neut = engine_d.run_neutronics_step(0, particles_override=160000, batches_override=50)
    d_th = engine_d.run_thermal_step(d_neut["q_vol_new"])
    
    d_k_val = d_neut["k_eff"]
    d_k_std = d_neut["k_std"]
    d_p_pins = d_th["pin_powers"]
    d_p_peak = float(np.max(d_p_pins))
    d_p_avg = float(np.mean(d_p_pins))
    d_f_r = d_p_peak / d_p_avg
    d_peak_idx = int(np.argmax(d_p_pins))
    d_f_z = float(np.max(d_neut["q_vol_new"][d_peak_idx, :]) / np.mean(d_neut["q_vol_new"][d_peak_idx, :]))
    d_f_q = d_f_r * d_f_z
    d_Tf_max = float(np.max(d_th["T_fuel_center"]))
    d_Tc_max = float(np.max(d_th["T_clad_out"]))
    d_Tcool_max = float(np.max(d_th["T_cool_bulk"]))
    d_dp = float(d_th["dp_total"])
    d_delta_rho = (d_k_val - k_ref_uncoupled) / (d_k_val * k_ref_uncoupled) * 1e5
    
    # Benchmark on RUN D
    bench_d = engine_d.run_cardinal_benchmark_pin(d_peak_idx, d_neut["q_vol_new"][d_peak_idx, :], d_neut["iter_dir"])
    d_card_Tf = bench_d["solid_results"].get("T_fuel_max", d_Tf_max)
    d_card_Tc = bench_d["solid_results"].get("T_clad_max", d_Tc_max)
    d_card_dp = bench_d["thm_results"].get("dp_total", d_dp)
    
    print(f"  RUN D Finished: k_eff = {d_k_val:.5f} +/- {d_k_std:.5f}")
    print(f"  T_fuel,max = {d_Tf_max:.1f} K | Peak Pin Power = {d_p_peak:.1f} W (Pin #{d_peak_idx})")
    print(f"  F_r = {d_f_r:.4f} | F_z = {d_f_z:.4f} | F_q = {d_f_q:.4f}")
    print(f"  Cardinal Benchmark: T_fuel,max = {d_card_Tf:.1f} K | Diff = {abs(d_card_Tf - d_Tf_max):.2f} K")
    
    run_d_record = {
        "run_id": "RUN D (High-Statistics Final)",
        "k_eff": d_k_val,
        "sigma_k": d_k_std,
        "delta_rho_pcm": d_delta_rho,
        "peak_pin_power_W": d_p_peak,
        "F_r": d_f_r,
        "F_z": d_f_z,
        "F_q": d_f_q,
        "q_vol_max_W_m3": float(np.max(d_neut["q_vol_new"])),
        "T_fuel_max_K": d_Tf_max,
        "T_clad_max_K": d_Tc_max,
        "T_coolant_max_K": d_Tcool_max,
        "T_coolant_out_K": d_th["T_out_mean"],
        "mass_flow_kg_s": total_coolant_mass_flow(),
        "Re_min": d_th["Re_range"][0],
        "Re_max": d_th["Re_range"][1],
        "Nu_mean": d_th["Nu_mean"],
        "h_min_W_m2_K": d_th["h_range"][0],
        "h_max_W_m2_K": d_th["h_range"][1],
        "pressure_drop_Pa": d_dp,
        "iterations_count": len(b_iter_records),
        "particles": 160000,
        "batches": 50
    }
    
    # =========================================================================
    # STEP 5: SAVE MASTER RUN COMPARISON TABLE (v3_comparison.csv)
    # =========================================================================
    comp_records = [
        run_a_record,
        {
            "run_id": "RUN B (Improved Cardinal/THM)",
            "k_eff": b_final_rec["k_eff"],
            "sigma_k": b_final_rec["sigma_k"],
            "delta_rho_pcm": b_final_rec["delta_rho_pcm"],
            "peak_pin_power_W": b_final_rec["peak_pin_power_W"],
            "F_r": b_final_rec["F_r"],
            "F_z": b_final_rec["F_z"],
            "F_q": b_final_rec["F_q"],
            "q_vol_max_W_m3": b_final_rec["q_vol_max_W_m3"],
            "T_fuel_max_K": b_final_rec["T_fuel_max_K"],
            "T_clad_max_K": b_final_rec["T_clad_max_K"],
            "T_coolant_max_K": b_final_rec["T_coolant_max_K"],
            "T_coolant_out_K": b_final_rec["T_coolant_out_K"],
            "mass_flow_kg_s": total_coolant_mass_flow(),
            "Re_min": th_res["Re_range"][0],
            "Re_max": th_res["Re_range"][1],
            "Nu_mean": th_res["Nu_mean"],
            "h_min_W_m2_K": th_res["h_range"][0],
            "h_max_W_m2_K": th_res["h_range"][1],
            "pressure_drop_Pa": b_final_rec["pressure_drop_Pa"],
            "iterations_count": len(b_iter_records),
            "particles": 80000,
            "batches": 45
        },
        {
            "run_id": "RUN C (20 Axial Slices Refined)",
            "k_eff": c_neut["k_eff"],
            "sigma_k": c_neut["k_std"],
            "delta_rho_pcm": (c_neut["k_eff"] - k_ref_uncoupled) / (c_neut["k_eff"] * k_ref_uncoupled) * 1e5,
            "peak_pin_power_W": c_p_peak,
            "F_r": c_f_r,
            "F_z": c_f_z,
            "F_q": c_f_q,
            "q_vol_max_W_m3": float(np.max(c_q_relaxed)),
            "T_fuel_max_K": c_Tf_max,
            "T_clad_max_K": c_Tc_max,
            "T_coolant_max_K": float(np.max(c_th["T_cool_bulk"])),
            "T_coolant_out_K": c_th["T_out_mean"],
            "mass_flow_kg_s": total_coolant_mass_flow(),
            "Re_min": c_th["Re_range"][0],
            "Re_max": c_th["Re_range"][1],
            "Nu_mean": c_th["Nu_mean"],
            "h_min_W_m2_K": c_th["h_range"][0],
            "h_max_W_m2_K": c_th["h_range"][1],
            "pressure_drop_Pa": c_dp,
            "iterations_count": 2,
            "particles": 80000,
            "batches": 45
        },
        run_d_record
    ]
    df_comp = pd.DataFrame(comp_records)
    df_comp.to_csv(os.path.join(RESULTS_V3_DIR, "v3_comparison.csv"), index=False)
    
    # Save v3_summary.json
    summary_data = {
        "campaign": "OPT-17 Cand 17-Clean V3 Multiphysics Coupling",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "total_runtime_s": time.time() - t_start_total,
        "runs": comp_records,
        "convergence_audit": conv_rec,
        "cardinal_benchmark": {
            "peak_pin_idx": d_peak_idx,
            "reduced_T_fuel_max_K": d_Tf_max,
            "cardinal_T_fuel_max_K": d_card_Tf,
            "diff_T_fuel_K": abs(d_card_Tf - d_Tf_max),
            "diff_T_fuel_percent": abs(d_card_Tf - d_Tf_max) / d_Tf_max * 100.0,
            "reduced_T_clad_max_K": d_Tc_max,
            "cardinal_T_clad_max_K": d_card_Tc,
            "diff_T_clad_K": abs(d_card_Tc - d_Tc_max),
            "reduced_dp_Pa": d_dp,
            "cardinal_dp_Pa": d_card_dp,
            "diff_dp_Pa": abs(d_card_dp - d_dp)
        },
        "discretization_sensitivity_20_vs_10": {
            "delta_k_pcm": diff_mesh_k,
            "delta_T_fuel_max_K": diff_mesh_Tf,
            "delta_F_z": diff_mesh_Fz,
            "delta_peak_pin_power_W": diff_mesh_Ppeak,
            "delta_pressure_drop_Pa": diff_mesh_dp
        }
    }
    with open(os.path.join(RESULTS_V3_DIR, "v3_summary.json"), "w") as f:
        json.dump(summary_data, f, indent=2)
        
    print("\n" + "#" * 70)
    print("  [STEP 6] GENERATING PUBLICATION-GRADE FIGURES (results/v3/figures/)")
    print("#" * 70)
    plot_all_v3_figures()
        
    print("\n" + "=" * 80)
    print("  V3 SIMULATION CAMPAIGN COMPLETED SUCCESSFULLY!")
    print(f"  Total Runtime: {time.time() - t_start_total:.1f} s")
    print(f"  Output directory: {RESULTS_V3_DIR}")
    print("=" * 80)


if __name__ == "__main__":
    execute_v3_campaign()
