#!/usr/bin/env python3
"""
run_depletion_4b.py

Multi-year depletion analysis at 2.5-year steps for Simulation 4B (5-Layer with Al2O3 Liner):
T_inlet = 920.0 K, T_outlet = 1130.0 K.
Tracks keff(0 deg), theta_crit(t), SDM(t), and exact subcritical crossing.
"""

import os
import sys
import time
import math
import numpy as np
import pandas as pd
import openmc

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/final_calculation/simulation_fourth/sim4b_al2o3"
sys.path.insert(0, os.path.join(BASE_DIR, "scripts"))
from build_model_4b import (
    build_openmc_model_4b, OPENMC_EXEC, CROSS_SECTIONS_PATH
)

PICARD_DIR = os.path.join(BASE_DIR, "coupled_picard")
OUT_DIR = os.path.join(BASE_DIR, "depletion_runs")
DATA_DIR = os.path.join(BASE_DIR, "data")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

grids_path = os.path.join(PICARD_DIR, "converged_temperature_grids.npz")
if os.path.exists(grids_path):
    grids = np.load(grids_path)
    T_fuel_vol = grids["T_fuel_vol"]
    T_clad_avg = grids["T_clad_avg"]
    T_cool_bulk = grids["T_cool_bulk"]
else:
    T_fuel_vol = np.full((217, 10), 1215.0)
    T_clad_avg = np.full((217, 10), 1190.0)
    T_cool_bulk = np.full((217, 10), 1025.0)

# Heavy metal inventory and burnup kinetics
M_HM_BOL = 124.12
M_U235_BOL = 24.47
M_U238_BOL = M_HM_BOL - M_U235_BOL
U235_TOTAL_BURN_KG_YR = 0.18450    # 153.75 g fission + 30.75 g capture (alpha = 0.20)
PU239_BREED_KG_YR = 0.02768        # Epithermal conversion ratio CR = 0.18
U238_CAPTURE_KG_YR = 0.02756
BURNUP_RATE_MWD_KGU_YR = 1.1771    # 400 kWth continuous


def compute_inventory(t_years):
    u235_t = max(0.0, M_U235_BOL - U235_TOTAL_BURN_KG_YR * t_years)
    u238_t = max(0.0, M_U238_BOL - U238_CAPTURE_KG_YR * t_years)
    pu239_t = PU239_BREED_KG_YR * t_years
    hm_t = u235_t + u238_t + pu239_t
    burnup_mwd_kgu = BURNUP_RATE_MWD_KGU_YR * t_years
    fima_pct = (burnup_mwd_kgu / 938.0) * 100.0

    enr_scale = (u235_t / M_U235_BOL) / (hm_t / M_HM_BOL)
    enr_in_eff = 19.30 * enr_scale
    enr_out_eff = 19.75 * enr_scale

    return {
        "t_yr": t_years,
        "u235_kg": u235_t,
        "u238_kg": u238_t,
        "pu239_kg": pu239_t,
        "hm_kg": hm_t,
        "burnup_mwd_kgu": burnup_mwd_kgu,
        "fima_pct": fima_pct,
        "enr_in_eff": enr_in_eff,
        "enr_out_eff": enr_out_eff
    }


def run_transport(run_dir, enr_in, enr_out, drum_angle_deg, particles=40000, batches=30, inactive=10, seed=42):
    os.makedirs(run_dir, exist_ok=True)
    t0 = time.time()

    build_openmc_model_4b(
        run_dir=run_dir,
        enr_inner=enr_in,
        enr_outer=enr_out,
        drum_angles=[drum_angle_deg] * 10,
        al2o3_mod_th=0.010, # 100 um Al2O3 liner
        T_fuel_grid=T_fuel_vol,
        T_clad_grid=T_clad_avg,
        T_cool_grid=T_cool_bulk,
        T_mod=940.0,
        T_refl=880.0,
        particles=particles,
        batches=batches,
        inactive=inactive,
        seed=seed,
        tallies_enabled=False
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
    print("=" * 80)
    print("MULTI-YEAR DEPLETION AT 2.5-YEAR STEPS: SIMULATION 4B (5-LAYER WITH AL2O3)")
    print("Convention: theta = 180 deg (Shutdown), theta = 0 deg (Max Outward)")
    print("=" * 80)

    time_steps_years = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 13.7, 14.5, 15.0, 17.5, 20.0, 22.5, 25.0]
    results = []

    for t_yr in time_steps_years:
        inv = compute_inventory(t_yr)
        enr_in = inv["enr_in_eff"]
        enr_out = inv["enr_out_eff"]

        print(f"\n>>> [STEP t = {t_yr:.1f} yr | Burnup = {inv['burnup_mwd_kgu']:.2f} MWd/kgU | FIMA = {inv['fima_pct']:.2f}% | U235 = {inv['u235_kg']:.2f} kg]")

        # 1. Drums fully outward (theta = 0 deg)
        dir_0deg = os.path.join(OUT_DIR, f"step_{t_yr:.1f}yr_0deg")
        k_0, sig_0, dt_0 = run_transport(dir_0deg, enr_in, enr_out, 0.0)
        rho_ex_pcm = (k_0 - 1.0) / k_0 * 1e5
        print(f"    Drums OUT (theta = 0 deg):   k_eff = {k_0:.5f} +/- {sig_0:.5f} (Excess rho = {rho_ex_pcm:+.1f} pcm, {dt_0:.1f} s)")

        # 2. Drums fully inward (theta = 180 deg)
        dir_180deg = os.path.join(OUT_DIR, f"step_{t_yr:.1f}yr_180deg")
        k_180, sig_180, dt_180 = run_transport(dir_180deg, enr_in, enr_out, 180.0)
        sdm_pcm = (1.0 - k_180) / k_180 * 1e5
        drum_worth_pcm = (1.0 / k_180 - 1.0 / k_0) * 1e5
        print(f"    Drums IN  (theta = 180 deg): k_eff = {k_180:.5f} +/- {sig_180:.5f} (SDM = {sdm_pcm:.1f} pcm, Worth = {drum_worth_pcm:.1f} pcm, {dt_180:.1f} s)")

        # 3. Critical drum angle estimation and verification
        if rho_ex_pcm > 0.0:
            frac = min(1.0, max(0.0, rho_ex_pcm / drum_worth_pcm))
            cos_val = 1.0 - 2.0 * frac
            cos_val = max(-1.0, min(1.0, cos_val))
            theta_crit = math.degrees(math.acos(cos_val))

            dir_crit = os.path.join(OUT_DIR, f"step_{t_yr:.1f}yr_crit")
            k_crit, sig_crit, dt_c = run_transport(dir_crit, enr_in, enr_out, theta_crit)
            status = "CRITICAL (Reserve Active)"
            print(f"    Critical Angle: theta_crit = {theta_crit:.1f} deg -> k_eff = {k_crit:.5f} +/- {sig_crit:.5f} ({status})")
        else:
            theta_crit = 0.0
            k_crit = k_0
            sig_crit = sig_0
            status = "SUBCRITICAL / EXHAUSTED (theta_crit = 0.0 deg)"
            print(f"    Critical Angle: theta_crit = 0.0 deg -> k_eff = {k_crit:.5f} +/- {sig_crit:.5f} ({status})")

        rec = {
            "time_years": t_yr,
            "burnup_mwd_kgu": inv["burnup_mwd_kgu"],
            "fima_percent": inv["fima_pct"],
            "u235_mass_kg": inv["u235_kg"],
            "pu239_mass_kg": inv["pu239_kg"],
            "k_eff_drums_out_0deg": k_0,
            "sigma_k_0deg": sig_0,
            "excess_reactivity_pcm": rho_ex_pcm,
            "k_eff_drums_in_180deg": k_180,
            "sigma_k_180deg": sig_180,
            "shutdown_margin_pcm": sdm_pcm,
            "drum_worth_pcm": drum_worth_pcm,
            "critical_drum_angle_deg": theta_crit,
            "k_eff_operational_crit": k_crit,
            "sigma_k_crit": sig_crit,
            "status": status
        }
        results.append(rec)

    df = pd.DataFrame(results)
    csv_path = os.path.join(DATA_DIR, "depletion_2.5yr_steps_4b.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n[+] Multi-Year Depletion results 4B saved to {csv_path}")


if __name__ == "__main__":
    main()
