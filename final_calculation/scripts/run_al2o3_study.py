#!/usr/bin/env python3
"""
run_al2o3_study.py

Comprehensive scientific study of Al2O3 (Alumina) layers in V7_CAND_BALANCED_ROBUST:
1. Literature review synthesis for thickness and properties
2. Neutronic impact in OpenMC across thickness range:
   - Moderator sleeve (between coolant and YH2): t = [0.0, 0.01, 0.10, 0.25, 0.50] mm
   - Cladding coating (between cladding and coolant): t = [0.10, 0.25] mm
3. Thermal and hydrogen depletion performance analysis:
   - Temperature distributions across pin cell
   - Hydrogen equilibrium dissociation pressure P_H2(T) in YH2
   - Peak clad temperature and creep margin comparison
"""

import os
import sys
import time
import math
import numpy as np
import pandas as pd
import openmc

BASE_DIR = "/home/aleri/OpenMC_projects/thermohydraulic/final_calculation_audited"
sys.path.insert(0, os.path.join(BASE_DIR, "scripts"))
from build_model_with_al2o3 import (
    build_al2o3_openmc_model, OPENMC_EXEC, CROSS_SECTIONS_PATH
)

OUT_DIR = os.path.join(BASE_DIR, "data/al2o3_study_runs")
os.makedirs(OUT_DIR, exist_ok=True)

# Thermophysical & Kinetic constants
Q_CORE_W = 400.0e3
N_PINS = 217
CORE_H_M = 0.4245
F_Q = 1.741             # Total power peaking factor
Q_PRIME_AVG = Q_CORE_W / (N_PINS * CORE_H_M)   # ~4342 W/m
Q_PRIME_PEAK = Q_PRIME_AVG * F_Q               # ~7560 W/m

# Radii [m]
R_FUEL_M = 0.00580
R_CLAD_OUT_M = 0.00628

# Thermal conductivities [W/(m K)]
K_FUEL = 21.0           # UN fuel at ~1200 K
K_CLAD = 48.0           # Nb-1Zr at ~1200 K
K_AL2O3_DENSE = 8.5     # Sintered alpha-Al2O3 at ~1000 K
K_AL2O3_TBC = 3.0       # Plasma-sprayed TBC / porous alumina at ~1000 K
H_COOL = 2200.0         # Turbulent helium convection [W/(m2 K)]
H_GAP = 10000.0         # Fuel-clad helium gap conductance [W/(m2 K)]

# Coolant operating temperatures [K]
T_COOL_IN = 950.0
T_COOL_OUT_MIXED = 1150.0
T_COOL_PEAK_CHAN_OUT = 1167.5
T_COOL_PEAK_ELEV = 1100.0   # Coolant temperature at peak heat flux elevation (z ~ 0.75 L)


def solve_pin_temperatures(al2o3_clad_th_mm=0.0, al2o3_mod_th_mm=0.0, k_al2o3=K_AL2O3_DENSE,
                           T_sink_mod=900.0):
    """
    Solves 1D radial conjugate heat transfer at the peak axial elevation:
    Returns dict of temperatures and thermal margins.
    """
    t_clad_al_m = al2o3_clad_th_mm * 1e-3
    t_mod_al_m = al2o3_mod_th_mm * 1e-3
    
    r_f = R_FUEL_M
    r_c = R_CLAD_OUT_M
    r_c_al = r_c + t_clad_al_m
    
    # 1. Convective drop to coolant from outer surface of pin (clad or clad+Al2O3)
    r_pin_out = r_c_al
    q_flux_pin_out = Q_PRIME_PEAK / (2.0 * math.pi * r_pin_out)
    delta_T_film = q_flux_pin_out / H_COOL
    T_surf = T_COOL_PEAK_ELEV + delta_T_film
    
    # 2. Temperature drop across inner Al2O3 layer (if present)
    if t_clad_al_m > 0.0:
        delta_T_clad_al = (Q_PRIME_PEAK / (2.0 * math.pi * k_al2o3)) * math.log(r_c_al / r_c)
    else:
        delta_T_clad_al = 0.0
        
    T_clad_out = T_surf + delta_T_clad_al
    
    # 3. Cladding conduction
    delta_T_clad = (Q_PRIME_PEAK / (2.0 * math.pi * K_CLAD)) * math.log(r_c / r_f)
    T_clad_in = T_clad_out + delta_T_clad
    
    # 4. Fuel-clad gap
    q_flux_fuel = Q_PRIME_PEAK / (2.0 * math.pi * r_f)
    delta_T_gap = q_flux_fuel / H_GAP
    T_fuel_surf = T_clad_in + delta_T_gap
    
    # 5. Fuel centerline conduction
    delta_T_fuel = Q_PRIME_PEAK / (4.0 * math.pi * K_FUEL)
    T_fuel_center = T_fuel_surf + delta_T_fuel
    
    # 6. Moderator temperature modeling
    # The moderator matrix receives heat from hot coolant and loses heat to reflector sink
    # Convective resistance at channel outer wall:
    r_ch_in = r_c_al + 0.00084  # coolant gap = 0.84 mm
    r_mod_inner = r_ch_in + t_mod_al_m
    
    R_conv_mod = 1.0 / (2.0 * math.pi * r_ch_in * H_COOL)
    R_ins_mod = (math.log(r_mod_inner / r_ch_in) / (2.0 * math.pi * k_al2o3)) if t_mod_al_m > 0 else 0.0
    R_sink_mod = 0.08   # Effective radial conduction / bypass resistance to reflector sink [K m / W]
    
    # Heat flow from coolant into moderator:
    # Q'_mod = (T_cool_avg - T_sink) / (R_conv + R_ins + R_sink)
    T_cool_avg = 1050.0
    R_total_mod = R_conv_mod + R_ins_mod + R_sink_mod
    q_prime_mod_leak = max(0.0, (T_cool_avg - T_sink_mod) / R_total_mod)
    
    # Moderator operating temperature:
    T_mod_oper = T_sink_mod + q_prime_mod_leak * R_sink_mod
    # Add nuclear volumetric heating offset (~15 K):
    T_mod_oper += 15.0
    
    # Hydrogen equilibrium dissociation pressure:
    # ln(P_H2 [atm]) = 23.82 - 25600 / T_mod
    p_h2_atm = math.exp(23.82 - 25600.0 / T_mod_oper)
    p_h2_kpa = p_h2_atm * 101.325
    
    # Permeation Reduction Factor (PRF):
    # Pure YH2 without Al2O3: PRF = 1
    # With Al2O3 coating: PRF = 1000 to 5000 (literature benchmark)
    prf = 1.0 if al2o3_mod_th_mm == 0.0 else (500.0 if al2o3_mod_th_mm < 0.05 else 2500.0)
    relative_h_leakage_rate = (p_h2_atm / prf)
    
    creep_limit_margin = 1350.0 - T_clad_out
    
    return {
        "al2o3_clad_th_mm": al2o3_clad_th_mm,
        "al2o3_mod_th_mm": al2o3_mod_th_mm,
        "k_al2o3": k_al2o3,
        "T_cool_peak_elev": T_COOL_PEAK_ELEV,
        "delta_T_clad_al": delta_T_clad_al,
        "T_clad_peak_K": T_clad_out,
        "creep_limit_margin_K": creep_limit_margin,
        "T_fuel_peak_K": T_fuel_center,
        "T_mod_oper_K": T_mod_oper,
        "p_h2_dissoc_atm": p_h2_atm,
        "p_h2_dissoc_kPa": p_h2_kpa,
        "prf": prf,
        "relative_h_leakage": relative_h_leakage_rate
    }


def run_neutronic_case(case_name, al2o3_clad_th_mm=0.0, al2o3_mod_th_mm=0.0,
                       particles=40000, batches=30, inactive=10, seed=42):
    """Executes OpenMC eigenvalue calculation for given geometry."""
    run_dir = os.path.join(OUT_DIR, case_name)
    os.makedirs(run_dir, exist_ok=True)
    t0 = time.time()
    
    info = build_al2o3_openmc_model(
        run_dir=run_dir,
        al2o3_clad_th=al2o3_clad_th_mm * 0.1,   # mm to cm
        al2o3_mod_th=al2o3_mod_th_mm * 0.1,     # mm to cm
        particles=particles,
        batches=batches,
        inactive=inactive,
        seed=seed
    )
    
    os.environ["OPENMC_CROSS_SECTIONS"] = CROSS_SECTIONS_PATH
    os.environ["OMP_NUM_THREADS"] = "8"
    cwd_prev = os.getcwd()
    os.chdir(run_dir)
    try:
        openmc.run(openmc_exec=OPENMC_EXEC, output=False)
    finally:
        os.chdir(cwd_prev)
        
    sp_path = os.path.join(run_dir, f"statepoint.{batches}.h5")
    with openmc.StatePoint(sp_path) as sp:
        k_val = float(sp.keff.nominal_value)
        k_sig = float(sp.keff.std_dev)
        
    dt = time.time() - t0
    return k_val, k_sig, dt


def main():
    print("=" * 85)
    print("      SCIENTIFIC STUDY OF AL2O3 LAYERS (MODERATOR SLEEVE & CLAD COATING)")
    print("=" * 85)
    
    cases = [
        {"cid": "CASE_BASE", "name": "Baseline Reference (No Al2O3)", "clad_th": 0.0, "mod_th": 0.0, "cat": "Baseline"},
        {"cid": "CASE_MOD_10um", "name": "Al2O3 Mod Liner 10 um (Thin HPB)", "clad_th": 0.0, "mod_th": 0.010, "cat": "Mod_HPB"},
        {"cid": "CASE_MOD_100um", "name": "Al2O3 Mod Liner 100 um (Thin TBC)", "clad_th": 0.0, "mod_th": 0.100, "cat": "Mod_TBC"},
        {"cid": "CASE_MOD_250um", "name": "Al2O3 Mod Liner 250 um (Ref Sleeve)", "clad_th": 0.0, "mod_th": 0.250, "cat": "Mod_Sleeve"},
        {"cid": "CASE_MOD_500um", "name": "Al2O3 Mod Liner 500 um (Heavy Sleeve)", "clad_th": 0.0, "mod_th": 0.500, "cat": "Mod_Heavy"},
        {"cid": "CASE_CLAD_100um", "name": "Al2O3 Clad Coating 100 um", "clad_th": 0.100, "mod_th": 0.0, "cat": "Clad_Coating"},
        {"cid": "CASE_CLAD_250um", "name": "Al2O3 Clad Coating 250 um", "clad_th": 0.250, "mod_th": 0.0, "cat": "Clad_Coating"},
        {"cid": "CASE_DUAL_250um", "name": "Dual Al2O3 (Clad 250 um + Mod 250 um)", "clad_th": 0.250, "mod_th": 0.250, "cat": "Dual_Layer"},
    ]
    
    # 1. Thermal Analysis for all cases
    print("\n--- 1. Conjugate Thermal & Hydrogen Dissociation Analysis ---")
    thermal_results = {}
    for c in cases:
        t_res = solve_pin_temperatures(al2o3_clad_th_mm=c["clad_th"], al2o3_mod_th_mm=c["mod_th"])
        thermal_results[c["cid"]] = t_res
        print(f"[{c['cid']:16s}] T_clad = {t_res['T_clad_peak_K']:.1f} K (Margin: {t_res['creep_limit_margin_K']:+5.1f} K) | "
              f"T_fuel = {t_res['T_fuel_peak_K']:.1f} K | T_mod = {t_res['T_mod_oper_K']:.1f} K | "
              f"P_H2 = {t_res['p_h2_dissoc_kPa']:.3e} kPa")

    # 2. Neutronic Runs in OpenMC
    print("\n--- 2. OpenMC Continuous-Energy Monte Carlo Neutronic Runs ---")
    rows = []
    k_base = None
    for c in cases:
        print(f">>> Running OpenMC for {c['name']} (clad={c['clad_th']} mm, mod={c['mod_th']} mm) ...")
        k_val, k_sig, dt = run_neutronic_case(c["cid"], c["clad_th"], c["mod_th"])
        if c["cid"] == "CASE_BASE":
            k_base = k_val
            delta_rho_pcm = 0.0
        else:
            delta_rho_pcm = ((k_val - k_base) / (k_val * k_base)) * 1e5
            
        t_res = thermal_results[c["cid"]]
        print(f"    -> keff = {k_val:.5f} +/- {k_sig:.5f} | Delta rho = {delta_rho_pcm:+6.1f} pcm | Runtime = {dt:.1f} s")
        
        rows.append({
            "case_id": c["cid"],
            "description": c["name"],
            "category": c["cat"],
            "al2o3_clad_th_mm": c["clad_th"],
            "al2o3_mod_th_mm": c["mod_th"],
            "k_eff": k_val,
            "sigma_k": k_sig,
            "reactivity_delta_pcm": delta_rho_pcm,
            "T_clad_peak_K": t_res["T_clad_peak_K"],
            "creep_margin_K": t_res["creep_limit_margin_K"],
            "T_fuel_peak_K": t_res["T_fuel_peak_K"],
            "T_mod_operating_K": t_res["T_mod_oper_K"],
            "p_h2_dissoc_kPa": t_res["p_h2_dissoc_kPa"],
            "permeation_reduction_factor": t_res["prf"],
            "relative_hydrogen_loss_rate": t_res["relative_h_leakage"],
            "runtime_s": dt
        })
        
    df = pd.DataFrame(rows)
    csv_path = os.path.join(BASE_DIR, "data/al2o3_sensitivity_and_thermal_study.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n[+] Saved complete Al2O3 study to:\n    {csv_path}")


if __name__ == "__main__":
    main()
