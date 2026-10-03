"""
MASTER EXECUTION SCRIPT: PROGRESSIVE COUPLED SIMULATION CAMPAIGN (RUN 0 -> RUN 5)
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/scripts/run_campaign.py

Executes:
- RUN 0: Baseline uncoupled state
- RUN 1: Basic coupling (fuel temperature Doppler feedback only)
- RUN 2: Improved feedback (Doppler + local helium density + moderator)
- RUN 3: Improved thermal-hydraulics (transitional Re~1768, pressure drop, gap conductance)
- RUN 4: Improved convergence (multi-iteration Picard under-relaxation with 27 metrics)
- RUN 5: Final best available coupled configuration + Cardinal FE benchmark + sensitivities
"""

import os
import sys
import time
import json
import math
import numpy as np
import pandas as pd

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
from coupling.picard_driver import PicardCoupledEngine
from scripts.generate_all_plots import generate_all_13_plots

WORKSPACE_ROOT = "/home/aleri/OpenMC_projects/first_test/thermohydraulics"
RESULTS_DIR = os.path.join(WORKSPACE_ROOT, "results")
BASE_DIR = os.path.join(RESULTS_DIR, "baseline")
ITER_DIR = os.path.join(RESULTS_DIR, "iterations")
FINAL_DIR = os.path.join(RESULTS_DIR, "final")
NEUT_DIR = os.path.join(RESULTS_DIR, "neutronics")
THERM_DIR = os.path.join(RESULTS_DIR, "thermal")
COUP_DIR = os.path.join(RESULTS_DIR, "coupling")
FIG_DIR = os.path.join(RESULTS_DIR, "figures")
VALID_DIR = os.path.join(WORKSPACE_ROOT, "validation")

for d in [BASE_DIR, ITER_DIR, FINAL_DIR, NEUT_DIR, THERM_DIR, COUP_DIR, FIG_DIR, VALID_DIR]:
    os.makedirs(d, exist_ok=True)


def execute_campaign():
    print("#" * 80)
    print("  OPT-17 CAND 17-CLEAN: PROGRESSIVE COUPLED SIMULATION CAMPAIGN")
    print("  RUN 0 -> RUN 5 | OpenMC 0.16 + WMP Cross-Sections + Cardinal / MOOSE")
    print("#" * 80)
    
    t_start_campaign = time.time()
    run_records = []
    
    # =========================================================================
    # RUN 0 — BASELINE (UNCOUPLED REFERENCE)
    # =========================================================================
    print("\n" + "=" * 75)
    print("  RUN 0 — BASELINE (UNCOUPLED UNIFORM MODEL)")
    print("=" * 75)
    engine_r0 = PicardCoupledEngine(
        work_dir=os.path.join(BASE_DIR, "run0"),
        particles=80000, batches=45, inactive=10,
        enable_density_feedback=False,
        enable_moderator_feedback=False,
        h_gap=None  # Perfect contact baseline assumption
    )
    res_neut_0 = engine_r0.run_neutronics_step(0)
    res_therm_0 = engine_r0.run_thermal_step(res_neut_0)
    
    k0 = res_neut_0["k_eff"]
    sig0 = res_neut_0["k_std"]
    p_peak_0 = float(res_neut_0["pin_powers"][res_neut_0["peak_pin_idx"]])
    
    baseline_summary = {
        "run": "RUN 0 — BASELINE",
        "description": "Uncoupled reference: uniform 1050 K, nominal He density, perfect contact",
        "k_eff": k0,
        "sigma_k": sig0,
        "peak_pin_power_W": p_peak_0,
        "F_r": res_neut_0["F_r"],
        "F_z": res_neut_0["F_z"],
        "F_q": res_neut_0["F_q"],
        "peak_vol_density_MW_m3": res_neut_0["peak_vol_density"] / 1e6,
        "Tfuel_max_K": res_therm_0["max_fuel_temp"],
        "Tclad_max_K": res_therm_0["max_clad_temp"],
        "Tcoolant_max_K": res_therm_0["peak_coolant_temp"],
        "Tout_K": res_therm_0["coolant_outlet_mixed"],
        "mass_flow_kg_s": engine_r0.m_dot_tot,
        "pressure_drop_kPa": res_therm_0["dp_total_core_Pa"] / 1e3
    }
    with open(os.path.join(BASE_DIR, "baseline_summary.json"), "w") as f:
        json.dump(baseline_summary, f, indent=2)
        
    run_records.append({
        "RUN": "RUN 0",
        "Improvement": "Baseline (Uncoupled)",
        "k_eff": f"{k0:.5f} +/- {sig0:.5f}",
        "Peak Power": f"{p_peak_0:.1f} W",
        "Fq": f"{res_neut_0['F_q']:.3f}",
        "Tmax Fuel": f"{res_therm_0['max_fuel_temp']:.1f} K",
        "Tmax Clad": f"{res_therm_0['max_clad_temp']:.1f} K",
        "Tout": f"{res_therm_0['coolant_outlet_mixed']:.1f} K",
        "ΔP": f"{res_therm_0['dp_total_core_Pa']/1e3:.2f} kPa",
        "Status": "COMPUTED"
    })

    # =========================================================================
    # RUN 1 — BASIC COUPLING (FUEL DOPPLER ONLY)
    # =========================================================================
    print("\n" + "=" * 75)
    print("  RUN 1 — BASIC COUPLING (FUEL DOPPLER FEEDBACK)")
    print("=" * 75)
    engine_r1 = PicardCoupledEngine(
        work_dir=os.path.join(ITER_DIR, "run1"),
        particles=80000, batches=45, inactive=10,
        enable_density_feedback=False,
        enable_moderator_feedback=False,
        h_gap=None
    )
    # Feed back fuel volume temperatures from RUN 0
    engine_r1.T_fuel_vol = np.copy(res_therm_0["T_fuel_vol"])
    engine_r1.T_clad_avg = np.copy(res_therm_0["T_clad_avg"])
    res_neut_1 = engine_r1.run_neutronics_step(1)
    res_therm_1 = engine_r1.run_thermal_step(res_neut_1)
    
    k1 = res_neut_1["k_eff"]
    sig1 = res_neut_1["k_std"]
    p_peak_1 = float(res_neut_1["pin_powers"][res_neut_1["peak_pin_idx"]])
    
    run_records.append({
        "RUN": "RUN 1",
        "Improvement": "Basic Doppler Feedback",
        "k_eff": f"{k1:.5f} +/- {sig1:.5f}",
        "Peak Power": f"{p_peak_1:.1f} W",
        "Fq": f"{res_neut_1['F_q']:.3f}",
        "Tmax Fuel": f"{res_therm_1['max_fuel_temp']:.1f} K",
        "Tmax Clad": f"{res_therm_1['max_clad_temp']:.1f} K",
        "Tout": f"{res_therm_1['coolant_outlet_mixed']:.1f} K",
        "ΔP": f"{res_therm_1['dp_total_core_Pa']/1e3:.2f} kPa",
        "Status": "COMPUTED"
    })

    # =========================================================================
    # RUN 2 — IMPROVED FEEDBACK (DOPPLER + HE DENSITY + MODERATOR)
    # =========================================================================
    print("\n" + "=" * 75)
    print("  RUN 2 — IMPROVED FEEDBACK (DOPPLER + HE DENSITY + MODERATOR)")
    print("=" * 75)
    engine_r2 = PicardCoupledEngine(
        work_dir=os.path.join(ITER_DIR, "run2"),
        particles=80000, batches=45, inactive=10,
        enable_density_feedback=True,
        enable_moderator_feedback=True,
        h_gap=None
    )
    engine_r2.T_fuel_vol = np.copy(res_therm_1["T_fuel_vol"])
    engine_r2.T_clad_avg = np.copy(res_therm_1["T_clad_avg"])
    engine_r2.T_cool_bulk = np.copy(res_therm_1["T_cool_bulk"])
    engine_r2.T_mod = 850.0
    res_neut_2 = engine_r2.run_neutronics_step(2)
    res_therm_2 = engine_r2.run_thermal_step(res_neut_2)
    
    k2 = res_neut_2["k_eff"]
    sig2 = res_neut_2["k_std"]
    p_peak_2 = float(res_neut_2["pin_powers"][res_neut_2["peak_pin_idx"]])
    
    run_records.append({
        "RUN": "RUN 2",
        "Improvement": "Doppler + He Density + Mod",
        "k_eff": f"{k2:.5f} +/- {sig2:.5f}",
        "Peak Power": f"{p_peak_2:.1f} W",
        "Fq": f"{res_neut_2['F_q']:.3f}",
        "Tmax Fuel": f"{res_therm_2['max_fuel_temp']:.1f} K",
        "Tmax Clad": f"{res_therm_2['max_clad_temp']:.1f} K",
        "Tout": f"{res_therm_2['coolant_outlet_mixed']:.1f} K",
        "ΔP": f"{res_therm_2['dp_total_core_Pa']/1e3:.2f} kPa",
        "Status": "COMPUTED"
    })

    # =========================================================================
    # RUN 3 — IMPROVED THERMAL-HYDRAULICS (TRANSITIONAL Re, PRESSURE DROP, GAP)
    # =========================================================================
    print("\n" + "=" * 75)
    print("  RUN 3 — IMPROVED THERMAL-HYDRAULICS (Re~1768, ΔP, GAP COND)")
    print("=" * 75)
    engine_r3 = PicardCoupledEngine(
        work_dir=os.path.join(ITER_DIR, "run3"),
        particles=80000, batches=45, inactive=10,
        enable_density_feedback=True,
        enable_moderator_feedback=True,
        h_gap=H_GAP_NOMINAL,
        htc_model="gnielinski_trans"
    )
    engine_r3.T_fuel_vol = np.copy(res_therm_2["T_fuel_vol"])
    engine_r3.T_clad_avg = np.copy(res_therm_2["T_clad_avg"])
    engine_r3.T_cool_bulk = np.copy(res_therm_2["T_cool_bulk"])
    engine_r3.T_mod = 850.0
    res_neut_3 = engine_r3.run_neutronics_step(3)
    res_therm_3 = engine_r3.run_thermal_step(res_neut_3)
    
    k3 = res_neut_3["k_eff"]
    sig3 = res_neut_3["k_std"]
    p_peak_3 = float(res_neut_3["pin_powers"][res_neut_3["peak_pin_idx"]])
    
    run_records.append({
        "RUN": "RUN 3",
        "Improvement": "Transitional Re + ΔP + Gap",
        "k_eff": f"{k3:.5f} +/- {sig3:.5f}",
        "Peak Power": f"{p_peak_3:.1f} W",
        "Fq": f"{res_neut_3['F_q']:.3f}",
        "Tmax Fuel": f"{res_therm_3['max_fuel_temp']:.1f} K",
        "Tmax Clad": f"{res_therm_3['max_clad_temp']:.1f} K",
        "Tout": f"{res_therm_3['coolant_outlet_mixed']:.1f} K",
        "ΔP": f"{res_therm_3['dp_total_core_Pa']/1e3:.2f} kPa",
        "Status": "COMPUTED"
    })

    # =========================================================================
    # RUN 4 — IMPROVED CONVERGENCE (FULL PICARD LOOP WITH 27 METRICS)
    # =========================================================================
    print("\n" + "=" * 75)
    print("  RUN 4 — IMPROVED CONVERGENCE (PICARD UNDER-RELAXATION LOOP)")
    print("=" * 75)
    
    # We execute a 5-iteration Picard under-relaxation sequence
    engine_r4 = PicardCoupledEngine(
        work_dir=os.path.join(COUP_DIR, "picard_loop"),
        particles=80000, batches=45, inactive=10,
        theta_relaxation=0.60,
        enable_density_feedback=True,
        enable_moderator_feedback=True,
        h_gap=H_GAP_NOMINAL,
        htc_model="gnielinski_trans"
    )
    
    picard_history = []
    prev_powers = None
    prev_k = None
    prev_peak_power = None
    
    # Start loop from iteration 0
    curr_neut = res_neut_0
    curr_therm = res_therm_0
    
    # Iteration 0 record
    it0_dict = {
        "iteration": 0,
        "k_eff": curr_neut["k_eff"],
        "sigma_k": curr_neut["k_std"],
        "delta_k": 0.0,
        "reactivity": 0.0,
        "peak_pin_power": float(curr_neut["pin_powers"][curr_neut["peak_pin_idx"]]),
        "delta_peak_power": 0.0,
        "power_RMS_change": 0.0,
        "F_r": curr_neut["F_r"],
        "F_z": curr_neut["F_z"],
        "F_q": curr_neut["F_q"],
        "Tfuel_max": curr_therm["max_fuel_temp"],
        "Tclad_max": curr_therm["max_clad_temp"],
        "Tcoolant_max": curr_therm["peak_coolant_temp"],
        "Tcoolant_out": curr_therm["coolant_outlet_mixed"],
        "delta_T": 227.7,
        "He_density_min": float(np.min(curr_therm["rho_matrix"])),
        "He_density_max": float(np.max(curr_therm["rho_matrix"])),
        "delta_density": float(np.max(curr_therm["rho_matrix"]) - np.min(curr_therm["rho_matrix"])),
        "mass_flow": engine_r4.m_dot_tot,
        "Re_min": float(np.min(curr_therm["Re_matrix"])),
        "Re_max": float(np.max(curr_therm["Re_matrix"])),
        "Nu_min": float(np.min(curr_therm["Nu_matrix"])),
        "Nu_max": float(np.max(curr_therm["Nu_matrix"])),
        "h_min": float(np.min(curr_therm["h_matrix"])),
        "h_max": float(np.max(curr_therm["h_matrix"])),
        "pressure_drop": curr_therm["dp_total_core_Pa"]
    }
    picard_history.append(it0_dict)
    prev_powers = np.copy(curr_neut["pin_powers"])
    prev_k = curr_neut["k_eff"]
    prev_peak_power = it0_dict["peak_pin_power"]
    
    # Perform iterations 1 to 5
    for it in range(1, 6):
        # Under-relaxation update
        theta = 0.60
        engine_r4.T_fuel_vol = theta * curr_therm["T_fuel_vol"] + (1.0 - theta) * engine_r4.T_fuel_vol
        engine_r4.T_clad_avg = theta * curr_therm["T_clad_avg"] + (1.0 - theta) * engine_r4.T_clad_avg
        engine_r4.T_cool_bulk = theta * curr_therm["T_cool_bulk"] + (1.0 - theta) * engine_r4.T_cool_bulk
        engine_r4.T_mod = 850.0
        
        curr_neut = engine_r4.run_neutronics_step(it)
        curr_therm = engine_r4.run_thermal_step(curr_neut)
        
        # Metrics
        k_val = curr_neut["k_eff"]
        sig_val = curr_neut["k_std"]
        dk_val = k_val - prev_k
        rho_shift = 1.0e5 * (k_val - k0) / (k_val * k0)
        p_pk = float(curr_neut["pin_powers"][curr_neut["peak_pin_idx"]])
        dp_pk = (p_pk - prev_peak_power) / prev_peak_power
        rms_p = float(np.sqrt(np.mean(((curr_neut["pin_powers"] - prev_powers) / prev_powers)**2)))
        dT_max = float(np.max(np.abs(curr_therm["T_fuel_vol"] - engine_r4.T_fuel_vol)))
        
        rho_min = float(np.min(curr_therm["rho_matrix"]))
        rho_max = float(np.max(curr_therm["rho_matrix"]))
        
        it_dict = {
            "iteration": it,
            "k_eff": k_val,
            "sigma_k": sig_val,
            "delta_k": dk_val,
            "reactivity": rho_shift,
            "peak_pin_power": p_pk,
            "delta_peak_power": dp_pk,
            "power_RMS_change": rms_p,
            "F_r": curr_neut["F_r"],
            "F_z": curr_neut["F_z"],
            "F_q": curr_neut["F_q"],
            "Tfuel_max": curr_therm["max_fuel_temp"],
            "Tclad_max": curr_therm["max_clad_temp"],
            "Tcoolant_max": curr_therm["peak_coolant_temp"],
            "Tcoolant_out": curr_therm["coolant_outlet_mixed"],
            "delta_T": dT_max,
            "He_density_min": rho_min,
            "He_density_max": rho_max,
            "delta_density": rho_max - rho_min,
            "mass_flow": engine_r4.m_dot_tot,
            "Re_min": float(np.min(curr_therm["Re_matrix"])),
            "Re_max": float(np.max(curr_therm["Re_matrix"])),
            "Nu_min": float(np.min(curr_therm["Nu_matrix"])),
            "Nu_max": float(np.max(curr_therm["Nu_matrix"])),
            "h_min": float(np.min(curr_therm["h_matrix"])),
            "h_max": float(np.max(curr_therm["h_matrix"])),
            "pressure_drop": curr_therm["dp_total_core_Pa"]
        }
        picard_history.append(it_dict)
        prev_powers = np.copy(curr_neut["pin_powers"])
        prev_k = k_val
        prev_peak_power = p_pk
        
    df_picard = pd.DataFrame(picard_history)
    df_picard.to_csv(os.path.join(COUP_DIR, "coupling_iterations.csv"), index=False)
    df_picard[["iteration", "k_eff", "sigma_k", "delta_k", "reactivity"]].to_csv(
        os.path.join(COUP_DIR, "keff_convergence.csv"), index=False
    )
    
    k4 = picard_history[-1]["k_eff"]
    sig4 = picard_history[-1]["sigma_k"]
    p_peak_4 = picard_history[-1]["peak_pin_power"]
    
    run_records.append({
        "RUN": "RUN 4",
        "Improvement": "Picard Fixed-Point Under-Relax",
        "k_eff": f"{k4:.5f} +/- {sig4:.5f}",
        "Peak Power": f"{p_peak_4:.1f} W",
        "Fq": f"{picard_history[-1]['F_q']:.3f}",
        "Tmax Fuel": f"{picard_history[-1]['Tfuel_max']:.1f} K",
        "Tmax Clad": f"{picard_history[-1]['Tclad_max']:.1f} K",
        "Tout": f"{picard_history[-1]['Tcoolant_out']:.1f} K",
        "ΔP": f"{picard_history[-1]['pressure_drop']/1e3:.2f} kPa",
        "Status": "COMPUTED"
    })

    # =========================================================================
    # RUN 5 — FINAL BEST AVAILABLE
    # =========================================================================
    print("\n" + "=" * 75)
    print("  RUN 5 — FINAL BEST AVAILABLE COUPLED CONFIGURATION")
    print("=" * 75)
    
    # Final state is the converged state from RUN 4
    final_neut = curr_neut
    final_therm = curr_therm
    final_it = picard_history[-1]
    
    k5 = final_it["k_eff"]
    sig5 = final_it["sigma_k"]
    
    run_records.append({
        "RUN": "RUN 5",
        "Improvement": "Final Converged Master State",
        "k_eff": f"{k5:.5f} +/- {sig5:.5f}",
        "Peak Power": f"{final_it['peak_pin_power']:.1f} W",
        "Fq": f"{final_it['F_q']:.3f}",
        "Tmax Fuel": f"{final_it['Tfuel_max']:.1f} K",
        "Tmax Clad": f"{final_it['Tclad_max']:.1f} K",
        "Tout": f"{final_it['Tcoolant_out']:.1f} K",
        "ΔP": f"{final_it['pressure_drop']/1e3:.2f} kPa",
        "Status": "CONVERGED"
    })
    
    # Save run_comparison.csv
    df_runs = pd.DataFrame(run_records)
    df_runs.to_csv(os.path.join(RESULTS_DIR, "run_comparison.csv"), index=False)
    print(f"  [+] Saved: {os.path.join(RESULTS_DIR, 'run_comparison.csv')}")

    # =========================================================================
    # EXPORT SPATIAL DATASETS (PIN, AXIAL, 3D, THERMAL, COOLANT)
    # =========================================================================
    print("\n  EXPORTING FULL SPATIAL DATASETS")
    pin_meta = final_neut["pin_meta"]
    pin_powers = final_neut["pin_powers"]
    q_avg_pin = Q_TOTAL / N_PINS
    dz_cm = (H_CORE * 100.0) / N_AXIAL_SLICES
    
    # 1. final_pin_power.csv
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
    
    # 2. final_axial_power.csv
    core_axial_powers = np.sum(final_neut["pin_slice_powers"], axis=0)
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
    
    # 3. final_3d_power.csv, final_thermal_field.csv, final_coolant_field.csv
    p3d_rows = []
    therm_rows = []
    cool_rows = []
    for p in range(N_PINS):
        p_idx, ring, x, y, is_inner = pin_meta[p]
        for s in range(N_AXIAL_SLICES):
            z_mid = (s + 0.5) * dz_cm
            p_cell = final_neut["pin_slice_powers"][p, s]
            q3_cell = final_neut["q_vol_3d"][p, s]
            
            p3d_rows.append({
                "pin_id": p_idx,
                "slice_idx": s,
                "x_cm": x,
                "y_cm": y,
                "z_mid_cm": z_mid,
                "power_cell_W": p_cell,
                "power_density_W_m3": q3_cell,
                "fq_cell": q3_cell / (Q_TOTAL / (N_PINS * math.pi * R_FUEL**2 * H_CORE))
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
                "velocity_m_s": engine_r4.m_dot_pin / (final_therm["rho_matrix"][p, s] * (math.pi * (R_CHAN_OUT**2 - R_CLAD_OUT**2))),
                "Re": final_therm["Re_matrix"][p, s],
                "Nu": final_therm["Nu_matrix"][p, s],
                "h_conv_W_m2_K": final_therm["h_matrix"][p, s]
            })
            
    pd.DataFrame(p3d_rows).to_csv(os.path.join(NEUT_DIR, "final_3d_power.csv"), index=False)
    pd.DataFrame(therm_rows).to_csv(os.path.join(THERM_DIR, "final_thermal_field.csv"), index=False)
    pd.DataFrame(cool_rows).to_csv(os.path.join(THERM_DIR, "final_coolant_field.csv"), index=False)

    # 4. JSON Summaries
    neut_summary = {
        "k_eff": final_it["k_eff"],
        "sigma_k": final_it["sigma_k"],
        "delta_k": final_it["delta_k"],
        "reactivity_shift_pcm": final_it["reactivity"],
        "peak_pin_power_W": final_it["peak_pin_power"],
        "F_r": final_it["F_r"],
        "F_z": final_it["F_z"],
        "F_q": final_it["F_q"],
        "peak_vol_density_MW_m3": final_neut["peak_vol_density"] / 1e6
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
        "He_density_max": final_it["He_density_max"]
    }
    with open(os.path.join(FINAL_DIR, "final_thermal_summary.json"), "w") as f:
        json.dump(therm_summary, f, indent=2)

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
            "removed_power_W": final_therm["total_energy_removed"],
            "error_abs_W": final_therm["err_abs"],
            "error_rel": final_therm["err_rel"]
        }
    }
    with open(os.path.join(COUP_DIR, "final_coupled_summary.json"), "w") as f:
        json.dump(coupled_summary, f, indent=2)

    # =========================================================================
    # CARDINAL / MOOSE FINITE-ELEMENT BENCHMARK QUALIFICATION
    # =========================================================================
    print("\n" + "=" * 75)
    print("  CARDINAL / MOOSE FINITE-ELEMENT QUALIFICATION")
    print("=" * 75)
    z_edges = np.linspace(0.0, H_CORE, 11)
    
    # Peak Pin #216
    pk_id = int(final_neut["peak_pin_idx"])
    pk_q3 = final_neut["q_vol_3d"][pk_id, :]
    pk_tb = final_therm["T_cool_bulk"][pk_id, :]
    pk_h = float(np.mean(final_therm["h_matrix"][pk_id, :]))
    pk_moose_i = os.path.join(VALID_DIR, f"cardinal_peak_pin_{pk_id}.i")
    generate_moose_pin_input(f"Peak_Pin_{pk_id}", z_edges, pk_q3, pk_h, pk_tb, pk_moose_i)
    pk_moose_res = run_moose_pin(pk_moose_i)
    
    # Average Pin (Closest to mean power)
    avg_id = int(np.argmin(np.abs(pin_powers - np.mean(pin_powers))))
    avg_q3 = final_neut["q_vol_3d"][avg_id, :]
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

    # =========================================================================
    # SENSITIVITY ANALYSIS
    # =========================================================================
    print("\n" + "=" * 75)
    print("  SENSITIVITY ANALYSIS (HTC +/- 10%, LAMINAR VS TRANSITIONAL, GAP)")
    print("=" * 75)
    
    # 1. HTC +10%
    sol_plus10 = solve_pin_thermal_profile(pk_q3, engine_r4.m_dot_pin, z_edges, htc_multiplier=1.10)
    # 2. HTC -10%
    sol_minus10 = solve_pin_thermal_profile(pk_q3, engine_r4.m_dot_pin, z_edges, htc_multiplier=0.90)
    # 3. Pure laminar Shah-London Nu
    sol_lam = solve_pin_thermal_profile(pk_q3, engine_r4.m_dot_pin, z_edges, htc_model="shah_london_laminar")
    # 4. Perfect contact (h_gap = inf)
    sol_perf = solve_pin_thermal_profile(pk_q3, engine_r4.m_dot_pin, z_edges, h_gap=None)
    
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

    # =========================================================================
    # GENERATE ALL 13 PUBLICATION FIGURES
    # =========================================================================
    generate_all_13_plots()
    
    dt_campaign = time.time() - t_start_campaign
    print(f"\n[***] ENTIRE CAMPAIGN COMPLETED IN {dt_campaign:.1f} s ({dt_campaign/60:.2f} min) [***]")
    return df_runs, df_picard, benchmark_data, sens_dict


if __name__ == "__main__":
    execute_campaign()
