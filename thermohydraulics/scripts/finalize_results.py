"""
FINALIZATION AND QUALIFICATION SCRIPT FOR OPT-17 CAMPAIGN
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/scripts/finalize_results.py

Processes the converged iteration 5 statepoint and thermal results:
1. Exports full spatial datasets (pin, axial, 3D power, thermal, coolant fields).
2. Performs Cardinal / MOOSE finite-element qualification for Peak Pin and Average Pin.
3. Performs sensitivity analyses (+/-10% HTC, pure laminar, gap conductance).
4. Generates all 13 publication figures at 300 DPI.
"""

import os
import sys
import json
import math
import numpy as np
import pandas as pd
import openmc

sys.path.insert(0, "/home/aleri/OpenMC_projects/first_test/thermohydraulics")
from models.th_properties import (
    Q_TOTAL, N_PINS, H_CORE, R_FUEL, R_CLAD_OUT, R_CHAN_OUT,
    T_INLET, DELTA_T_BULK, cp_helium, pin_coolant_mass_flow,
    total_coolant_mass_flow, rho_helium, H_GAP_NOMINAL,
    pressure_drop_channel
)
from models.moose_th_solver import (
    solve_pin_thermal_profile, generate_moose_pin_input, run_moose_pin
)
from inputs.build_opt17_openmc import (
    build_opt17_openmc_model, N_AXIAL_SLICES
)
from scripts.generate_all_plots import generate_all_13_plots

WORKSPACE_ROOT = "/home/aleri/OpenMC_projects/first_test/thermohydraulics"
RESULTS_DIR = os.path.join(WORKSPACE_ROOT, "results")
BASE_DIR = os.path.join(RESULTS_DIR, "baseline")
FINAL_DIR = os.path.join(RESULTS_DIR, "final")
NEUT_DIR = os.path.join(RESULTS_DIR, "neutronics")
THERM_DIR = os.path.join(RESULTS_DIR, "thermal")
COUP_DIR = os.path.join(RESULTS_DIR, "coupling")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
VALID_DIR = os.path.join(WORKSPACE_ROOT, "validation")

for d in [BASE_DIR, FINAL_DIR, NEUT_DIR, THERM_DIR, COUP_DIR, FIG_DIR, VALID_DIR]:
    os.makedirs(d, exist_ok=True)


def finalize_campaign():
    print("=" * 75)
    print("  FINALIZING OPT-17 CAMPAIGN DATASETS, BENCHMARKS, AND FIGURES")
    print("=" * 75)
    
    # 1. Load coupling_iterations.csv
    hist_csv = os.path.join(COUP_DIR, "coupling_iterations.csv")
    df_hist = pd.read_csv(hist_csv)
    final_it = df_hist.iloc[-1].to_dict()
    
    # 2. Extract final state from iter_5 statepoint
    iter5_dir = os.path.join(COUP_DIR, "picard_loop", "iter_5")
    sp_path = os.path.join(iter5_dir, "statepoint.45.h5")
    
    # Rebuild metadata to get cell maps
    model_info = build_opt17_openmc_model(
        run_dir=os.path.join(WORKSPACE_ROOT, "inputs", "temp_meta"),
        particles=100, batches=2, inactive=1
    )
    pin_meta = model_info["pin_meta"]
    cell_map = model_info["fuel_cells_map"]
    id_to_coords = {c.id: coords for coords, c in cell_map.items()}
    
    with openmc.StatePoint(sp_path) as sp:
        keff = sp.keff
        k_val = float(keff.nominal_value)
        k_std = float(keff.std_dev)
        tally = sp.get_tally(name="fuel_cell_heating")
        df_tally = tally.get_pandas_dataframe()
        
    df_heat = df_tally[df_tally["score"] == "heating"]
    raw_powers = np.zeros((N_PINS, N_AXIAL_SLICES))
    for _, row in df_heat.iterrows():
        cid = int(row["cell"])
        if cid in id_to_coords:
            p_idx, s_idx = id_to_coords[cid]
            raw_powers[p_idx, s_idx] = float(row["mean"])
            
    norm_factor = Q_TOTAL / np.sum(raw_powers)
    pin_slice_powers = raw_powers * norm_factor
    pin_powers = np.sum(pin_slice_powers, axis=1)
    
    v_fuel_slice = math.pi * (R_FUEL**2) * (H_CORE / N_AXIAL_SLICES)
    v_fuel_core = N_PINS * math.pi * (R_FUEL**2) * H_CORE
    q_vol_3d = pin_slice_powers / v_fuel_slice
    
    # Solve final thermal state with full transitional TH and gap conductance
    m_dot_pin = pin_coolant_mass_flow()
    z_edges = np.linspace(0.0, H_CORE, N_AXIAL_SLICES + 1)
    
    T_calc_fuel_vol = np.zeros((N_PINS, N_AXIAL_SLICES))
    T_calc_fuel_center = np.zeros((N_PINS, N_AXIAL_SLICES))
    T_calc_clad_out = np.zeros((N_PINS, N_AXIAL_SLICES))
    T_calc_clad_avg = np.zeros((N_PINS, N_AXIAL_SLICES))
    T_calc_cool_bulk = np.zeros((N_PINS, N_AXIAL_SLICES))
    h_matrix = np.zeros((N_PINS, N_AXIAL_SLICES))
    Re_matrix = np.zeros((N_PINS, N_AXIAL_SLICES))
    Nu_matrix = np.zeros((N_PINS, N_AXIAL_SLICES))
    rho_matrix = np.zeros((N_PINS, N_AXIAL_SLICES))
    dp_fric_arr = np.zeros(N_PINS)
    dp_accel_arr = np.zeros(N_PINS)
    
    for p in range(N_PINS):
        sol_p = solve_pin_thermal_profile(
            q_vol_3d[p, :], m_dot_pin, z_edges,
            h_gap=H_GAP_NOMINAL,
            htc_model="gnielinski_trans"
        )
        T_calc_fuel_vol[p, :] = sol_p["T_fuel_vol"]
        T_calc_fuel_center[p, :] = sol_p["T_fuel_center"]
        T_calc_clad_out[p, :] = sol_p["T_clad_out"]
        T_calc_clad_avg[p, :] = 0.5 * (sol_p["T_clad_out"] + sol_p["T_clad_in"])
        T_calc_cool_bulk[p, :] = sol_p["T_bulk"]
        h_matrix[p, :] = sol_p["h_conv"]
        Re_matrix[p, :] = sol_p["Re"]
        Nu_matrix[p, :] = sol_p["Nu"]
        rho_matrix[p, :] = sol_p["rho_cool"]
        dp_fric_arr[p] = sol_p["dp_fric_total_Pa"]
        dp_accel_arr[p] = sol_p["dp_accel_total_Pa"]
        
    final_therm = {
        "T_fuel_vol": T_calc_fuel_vol,
        "T_fuel_center": T_calc_fuel_center,
        "T_clad_out": T_calc_clad_out,
        "T_clad_avg": T_calc_clad_avg,
        "T_cool_bulk": T_calc_cool_bulk,
        "h_matrix": h_matrix,
        "Re_matrix": Re_matrix,
        "Nu_matrix": Nu_matrix,
        "rho_matrix": rho_matrix,
        "max_fuel_temp": float(np.max(T_calc_fuel_center)),
        "max_clad_temp": float(np.max(T_calc_clad_out)),
        "coolant_outlet_mixed": float(np.mean(T_calc_cool_bulk[:, -1])),
        "peak_coolant_temp": float(np.max(T_calc_cool_bulk[:, -1])),
        "dp_fric_core": float(np.mean(dp_fric_arr)),
        "dp_accel_core": float(np.mean(dp_accel_arr))
    }
    
    # 3. Export Spatial Datasets
    print("  [+] Exporting final_pin_power.csv...")
    q_avg_pin = Q_TOTAL / N_PINS
    dz_cm = (H_CORE * 100.0) / N_AXIAL_SLICES
    
    pin_rows = []
    for p in range(N_PINS):
        p_idx, ring, x, y, is_inner = pin_meta[p]
        p_val = pin_powers[p]
        pin_rows.append({
            "pin_id": p_idx,
            "ring": ring,
            "x_cm": x,
            "y_cm": y,
            "is_inner_zone": is_inner,
            "power_W": p_val,
            "radial_peaking_Fr": p_val / q_avg_pin,
            "T_coolant_exit_K": final_therm["T_cool_bulk"][p, -1],
            "T_fuel_max_center_K": np.max(final_therm["T_fuel_center"][p, :]),
            "T_clad_max_outer_K": np.max(final_therm["T_clad_out"][p, :])
        })
    pd.DataFrame(pin_rows).to_csv(os.path.join(NEUT_DIR, "final_pin_power.csv"), index=False)
    
    print("  [+] Exporting final_axial_power.csv...")
    core_axial_powers = np.sum(pin_slice_powers, axis=0)
    q_avg_slice = Q_TOTAL / N_AXIAL_SLICES
    axial_rows = []
    for s in range(N_AXIAL_SLICES):
        z_lo = s * dz_cm
        z_hi = (s + 1) * dz_cm
        z_mid = 0.5 * (z_lo + z_hi)
        axial_rows.append({
            "slice_idx": s,
            "z_start_cm": z_lo,
            "z_end_cm": z_hi,
            "z_mid_cm": z_mid,
            "slice_power_W": core_axial_powers[s],
            "axial_peaking_Fz": core_axial_powers[s] / q_avg_slice,
            "linear_power_W_cm": core_axial_powers[s] / (dz_cm * N_PINS),
            "coolant_bulk_mean_K": float(np.mean(final_therm["T_cool_bulk"][:, s])),
            "fuel_vol_mean_K": float(np.mean(final_therm["T_fuel_vol"][:, s])),
            "fuel_center_mean_K": float(np.mean(final_therm["T_fuel_center"][:, s]))
        })
    pd.DataFrame(axial_rows).to_csv(os.path.join(NEUT_DIR, "final_axial_power.csv"), index=False)
    
    print("  [+] Exporting final_3d_power.csv, final_thermal_field.csv, final_coolant_field.csv...")
    p3d_rows = []
    therm_rows = []
    cool_rows = []
    for p in range(N_PINS):
        p_idx, ring, x, y, is_inner = pin_meta[p]
        for s in range(N_AXIAL_SLICES):
            z_mid = (s + 0.5) * dz_cm
            p_cell = pin_slice_powers[p, s]
            q3_cell = q_vol_3d[p, s]
            
            p3d_rows.append({
                "pin_id": p_idx,
                "slice_idx": s,
                "x_cm": x,
                "y_cm": y,
                "z_mid_cm": z_mid,
                "power_cell_W": p_cell,
                "power_density_W_m3": q3_cell,
                "fq_cell": q3_cell / (Q_TOTAL / v_fuel_core)
            })
            
            therm_rows.append({
                "pin_id": p_idx,
                "slice_idx": s,
                "x_cm": x,
                "y_cm": y,
                "z_mid_cm": z_mid,
                "T_cool_bulk_K": final_therm["T_cool_bulk"][p, s],
                "T_clad_out_K": final_therm["T_clad_out"][p, s],
                "T_clad_avg_K": final_therm["T_clad_avg"][p, s],
                "T_fuel_vol_K": final_therm["T_fuel_vol"][p, s],
                "T_fuel_center_K": final_therm["T_fuel_center"][p, s],
                "h_conv_W_m2_K": final_therm["h_matrix"][p, s]
            })
            
            cool_rows.append({
                "pin_id": p_idx,
                "slice_idx": s,
                "z_mid_cm": z_mid,
                "T_bulk_K": final_therm["T_cool_bulk"][p, s],
                "density_kg_m3": final_therm["rho_matrix"][p, s],
                "velocity_m_s": m_dot_pin / (final_therm["rho_matrix"][p, s] * (math.pi * (R_CHAN_OUT**2 - R_CLAD_OUT**2))),
                "Re": final_therm["Re_matrix"][p, s],
                "Nu": final_therm["Nu_matrix"][p, s],
                "h_conv_W_m2_K": final_therm["h_matrix"][p, s]
            })
            
    pd.DataFrame(p3d_rows).to_csv(os.path.join(NEUT_DIR, "final_3d_power.csv"), index=False)
    pd.DataFrame(therm_rows).to_csv(os.path.join(THERM_DIR, "final_thermal_field.csv"), index=False)
    pd.DataFrame(cool_rows).to_csv(os.path.join(THERM_DIR, "final_coolant_field.csv"), index=False)

    # 4. JSON Summaries
    print("  [+] Exporting JSON summaries...")
    neut_summary = {
        "k_eff": final_it["k_eff"],
        "sigma_k": final_it["sigma_k"],
        "delta_k": final_it["delta_k"],
        "reactivity_shift_pcm": final_it["reactivity"],
        "peak_pin_power_W": final_it["peak_pin_power"],
        "F_r": final_it["F_r"],
        "F_z": final_it["F_z"],
        "F_q": final_it["F_q"],
        "peak_vol_density_MW_m3": float(np.max(q_vol_3d)) / 1e6
    }
    with open(os.path.join(FINAL_DIR, "final_neutronics_summary.json"), "w") as f:
        json.dump(neut_summary, f, indent=2)

    therm_summary = {
        "Tfuel_max_K": final_it["Tfuel_max"],
        "Tclad_max_K": final_it["Tclad_max"],
        "Tcoolant_max_K": final_it["Tcoolant_max"],
        "Tcoolant_out_K": final_it["Tcoolant_out"],
        "pressure_drop_kPa": final_it["pressure_drop"] / 1e3,
        "Re_min": final_it["Re_min"],
        "Re_max": final_it["Re_max"],
        "Nu_min": final_it["Nu_min"],
        "Nu_max": final_it["Nu_max"],
        "h_min": final_it["h_min"],
        "h_max": final_it["h_max"],
        "He_density_min": final_it["He_density_min"],
        "He_density_max": final_it["He_density_max"],
        "mass_flow_kg_s": total_coolant_mass_flow()
    }
    with open(os.path.join(FINAL_DIR, "final_thermal_summary.json"), "w") as f:
        json.dump(therm_summary, f, indent=2)

    with open(os.path.join(BASE_DIR, "baseline_summary.json")) as f:
        baseline_summary = json.load(f)

    coupled_summary = {
        "reactor": "Opt-17 Cand 17-Clean",
        "power_kW": Q_TOTAL / 1e3,
        "baseline_uncoupled": baseline_summary,
        "final_coupled": {
            "neutronics": neut_summary,
            "thermal": therm_summary
        },
        "energy_balance": {
            "generated_power_W": Q_TOTAL,
            "removed_power_W": float(np.sum(m_dot_pin * cp_helium() * (final_therm["T_cool_bulk"][:, -1] - T_INLET))),
            "error_abs_W": float(abs(np.sum(m_dot_pin * cp_helium() * (final_therm["T_cool_bulk"][:, -1] - T_INLET)) - Q_TOTAL)),
            "error_rel": float(abs(np.sum(m_dot_pin * cp_helium() * (final_therm["T_cool_bulk"][:, -1] - T_INLET)) - Q_TOTAL) / Q_TOTAL)
        }
    }
    with open(os.path.join(COUP_DIR, "final_coupled_summary.json"), "w") as f:
        json.dump(coupled_summary, f, indent=2)

    # 5. Cardinal / MOOSE Finite Element Benchmark
    print("\n  [+] Executing Cardinal / MOOSE FE Benchmarks on Peak and Average Pins...")
    pk_id = int(np.argmax(pin_powers))
    pk_q3 = q_vol_3d[pk_id, :]
    pk_tb = final_therm["T_cool_bulk"][pk_id, :]
    pk_h = float(np.mean(final_therm["h_matrix"][pk_id, :]))
    pk_moose_i = os.path.join(VALID_DIR, f"cardinal_peak_pin_{pk_id}.i")
    generate_moose_pin_input(f"Peak_Pin_{pk_id}", z_edges, pk_q3, pk_h, pk_tb, pk_moose_i)
    pk_moose_res = run_moose_pin(pk_moose_i)
    
    avg_id = int(np.argmin(np.abs(pin_powers - np.mean(pin_powers))))
    avg_q3 = q_vol_3d[avg_id, :]
    avg_tb = final_therm["T_cool_bulk"][avg_id, :]
    avg_h = float(np.mean(final_therm["h_matrix"][avg_id, :]))
    avg_moose_i = os.path.join(VALID_DIR, f"cardinal_avg_pin_{avg_id}.i")
    generate_moose_pin_input(f"Avg_Pin_{avg_id}", z_edges, avg_q3, avg_h, avg_tb, avg_moose_i)
    avg_moose_res = run_moose_pin(avg_moose_i)
    
    benchmark_data = {
        "peak_pin": {
            "pin_id": pk_id,
            "power_W": float(pin_powers[pk_id]),
            "model_Tfuel_max": float(np.max(final_therm["T_fuel_center"][pk_id, :])),
            "cardinal_Tfuel_max": float(pk_moose_res["T_fuel_max"]),
            "diff_Tfuel_max": float(abs(pk_moose_res["T_fuel_max"] - np.max(final_therm["T_fuel_center"][pk_id, :]))),
            "model_Tclad_max": float(np.max(final_therm["T_clad_out"][pk_id, :])),
            "cardinal_Tclad_max": float(pk_moose_res["T_clad_max"]),
            "diff_Tclad_max": float(abs(pk_moose_res["T_clad_max"] - np.max(final_therm["T_clad_out"][pk_id, :]))),
            "cardinal_Q_removed": float(pk_moose_res["Q_removed_boundary"])
        },
        "average_pin": {
            "pin_id": avg_id,
            "power_W": float(pin_powers[avg_id]),
            "model_Tfuel_max": float(np.max(final_therm["T_fuel_center"][avg_id, :])),
            "cardinal_Tfuel_max": float(avg_moose_res["T_fuel_max"]),
            "diff_Tfuel_max": float(abs(avg_moose_res["T_fuel_max"] - np.max(final_therm["T_fuel_center"][avg_id, :]))),
            "model_Tclad_max": float(np.max(final_therm["T_clad_out"][avg_id, :])),
            "cardinal_Tclad_max": float(avg_moose_res["T_clad_max"]),
            "diff_Tclad_max": float(abs(avg_moose_res["T_clad_max"] - np.max(final_therm["T_clad_out"][avg_id, :]))),
            "cardinal_Q_removed": float(avg_moose_res["Q_removed_boundary"])
        }
    }
    with open(os.path.join(VALID_DIR, "cardinal_benchmark_results.json"), "w") as f:
        json.dump(benchmark_data, f, indent=2)

    # 6. Sensitivity Analysis
    print("  [+] Performing Sensitivity Studies (+/-10% HTC, laminar vs transitional, gap conductance)...")
    sol_plus10 = solve_pin_thermal_profile(pk_q3, m_dot_pin, z_edges, htc_multiplier=1.10)
    sol_minus10 = solve_pin_thermal_profile(pk_q3, m_dot_pin, z_edges, htc_multiplier=0.90)
    sol_lam = solve_pin_thermal_profile(pk_q3, m_dot_pin, z_edges, htc_model="shah_london_laminar")
    sol_perf = solve_pin_thermal_profile(pk_q3, m_dot_pin, z_edges, h_gap=None)
    
    sens_dict = {
        "nominal_transitional_gap": {
            "Tfuel_max": float(np.max(final_therm["T_fuel_center"][pk_id, :])),
            "Tclad_max": float(np.max(final_therm["T_clad_out"][pk_id, :])),
            "h_mean": float(np.mean(final_therm["h_matrix"][pk_id, :]))
        },
        "htc_plus_10_percent": {
            "Tfuel_max": float(np.max(sol_plus10["T_fuel_center"])),
            "Tclad_max": float(np.max(sol_plus10["T_clad_out"])),
            "h_mean": float(np.mean(sol_plus10["h_conv"])),
            "delta_Tfuel": float(np.max(sol_plus10["T_fuel_center"]) - np.max(final_therm["T_fuel_center"][pk_id, :]))
        },
        "htc_minus_10_percent": {
            "Tfuel_max": float(np.max(sol_minus10["T_fuel_center"])),
            "Tclad_max": float(np.max(sol_minus10["T_clad_out"])),
            "h_mean": float(np.mean(sol_minus10["h_conv"])),
            "delta_Tfuel": float(np.max(sol_minus10["T_fuel_center"]) - np.max(final_therm["T_fuel_center"][pk_id, :]))
        },
        "shah_london_pure_laminar": {
            "Tfuel_max": float(np.max(sol_lam["T_fuel_center"])),
            "Tclad_max": float(np.max(sol_lam["T_clad_out"])),
            "h_mean": float(np.mean(sol_lam["h_conv"])),
            "delta_Tfuel": float(np.max(sol_lam["T_fuel_center"]) - np.max(final_therm["T_fuel_center"][pk_id, :]))
        },
        "perfect_fuel_clad_contact": {
            "Tfuel_max": float(np.max(sol_perf["T_fuel_center"])),
            "Tclad_max": float(np.max(sol_perf["T_clad_out"])),
            "delta_Tfuel": float(np.max(sol_perf["T_fuel_center"]) - np.max(final_therm["T_fuel_center"][pk_id, :]))
        }
    }
    with open(os.path.join(VALID_DIR, "sensitivity_results.json"), "w") as f:
        json.dump(sens_dict, f, indent=2)

    # 7. Generate All 13 Figures
    generate_all_13_plots()
    
    print("\n[***] ALL DATASETS, BENCHMARKS, SENSITIVITIES, AND 13 PLOTS FINALIZED! [***]")


if __name__ == "__main__":
    finalize_campaign()
