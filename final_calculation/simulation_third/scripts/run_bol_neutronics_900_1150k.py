#!/usr/bin/env python3
"""
run_bol_neutronics_900_1150k.py

High-statistics Beginning-of-Life (BOL) neutronics suite for Third Simulation (T_inlet = 900.0 K, T_outlet = 1150.0 K):
Evaluates:
  1. Operating state with all drums fully withdrawn (theta = 0.0 deg)
  2. Normal shutdown with all drums inserted (theta = 180.0 deg)
  3. Single stuck-drum safety margin (1 drum at 0 deg, 9 drums at 180 deg)
  4. Operating critical drum angle theta_crit where k_eff = 1.0000
Uses converged temperature fields from coupled Picard iteration.
"""

import os
import sys
import time
import math
import numpy as np
import pandas as pd
import openmc

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/final_calculation/simulation_third"
sys.path.insert(0, os.path.join(BASE_DIR, "scripts"))
from build_reference_model_900_1150k import (
    build_reference_openmc_model, OPENMC_EXEC, CROSS_SECTIONS_PATH
)

PICARD_DIR = os.path.join(BASE_DIR, "coupled_picard")
OUT_DIR = os.path.join(BASE_DIR, "bol_neutronics")
DATA_DIR = os.path.join(BASE_DIR, "data")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# Load converged temperature grids if available, else use defaults
grids_path = os.path.join(PICARD_DIR, "converged_temperature_grids.npz")
if os.path.exists(grids_path):
    grids = np.load(grids_path)
    T_fuel_vol = grids["T_fuel_vol"]
    T_clad_avg = grids["T_clad_avg"]
    T_cool_bulk = grids["T_cool_bulk"]
    print(f"[+] Loaded converged temperature grids from {grids_path}")
else:
    T_fuel_vol = np.full((217, 10), 1210.0)
    T_clad_avg = np.full((217, 10), 1195.0)
    T_cool_bulk = np.full((217, 10), 1025.0)
    print("[-] Using nominal 900-1150 K temperature grids")


def run_state(case_name, drum_angles, particles=50000, batches=40, inactive=10, seed=42):
    run_dir = os.path.join(OUT_DIR, case_name)
    os.makedirs(run_dir, exist_ok=True)
    t0 = time.time()

    build_reference_openmc_model(
        run_dir=run_dir,
        drum_angles=drum_angles,
        T_fuel_grid=T_fuel_vol,
        T_clad_grid=T_clad_avg,
        T_cool_grid=T_cool_bulk,
        T_mod=940.0,
        T_refl=875.0,
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
    print("RUNNING BOL HIGH-STATISTICS NEUTRONICS SUITE (T_inlet = 900.0 K, T_outlet = 1150.0 K)")
    print("=" * 80)

    # 1. Operating State: theta = 0.0 deg (All drums out)
    print("\n[1/4] Running Operating State (All Drums Out, theta = 0.0 deg)...")
    k_0, sig_0, dt_0 = run_state("bol_operating_0deg", [0.0] * 10)
    rho_ex_pcm = (k_0 - 1.0) / k_0 * 1e5
    sig_rho_ex = (sig_0 / (k_0**2)) * 1e5
    print(f"    k_eff(0 deg) = {k_0:.5f} +/- {sig_0:.5f} (Excess rho = {rho_ex_pcm:+.1f} +/- {sig_rho_ex:.1f} pcm, {dt_0:.1f} s)")

    # 2. Normal Shutdown: theta = 180.0 deg (All drums in)
    print("\n[2/4] Running Normal Shutdown (All Drums In, theta = 180.0 deg)...")
    k_180, sig_180, dt_180 = run_state("bol_shutdown_180deg", [180.0] * 10)
    sdm_all_pcm = (1.0 - k_180) / k_180 * 1e5
    sig_sdm_all = (sig_180 / (k_180**2)) * 1e5
    drum_worth_pcm = (1.0 / k_180 - 1.0 / k_0) * 1e5
    print(f"    k_eff(180 deg) = {k_180:.5f} +/- {sig_180:.5f} (SDM = {sdm_all_pcm:.1f} +/- {sig_sdm_all:.1f} pcm, {dt_180:.1f} s)")
    print(f"    Total Drum Worth = {drum_worth_pcm:.1f} pcm")

    # 3. Stuck Drum: 1 drum at 0 deg, 9 drums at 180 deg
    print("\n[3/4] Running Single Stuck Drum (1 Drum OUT, 9 Drums IN)...")
    stuck_angles = [0.0] + [180.0] * 9
    k_stuck, sig_stuck, dt_stuck = run_state("bol_stuck_drum", stuck_angles)
    sdm_stuck_pcm = (1.0 - k_stuck) / k_stuck * 1e5
    sig_sdm_stuck = (sig_stuck / (k_stuck**2)) * 1e5
    print(f"    k_eff(stuck) = {k_stuck:.5f} +/- {sig_stuck:.5f} (SDM_stuck = {sdm_stuck_pcm:.1f} +/- {sig_sdm_stuck:.1f} pcm, {dt_stuck:.1f} s)")

    # 4. Critical Drum Angle Search
    print("\n[4/4] Evaluating BOL Critical Drum Angle...")
    frac = rho_ex_pcm / drum_worth_pcm
    theta_crit_est = math.degrees(math.acos(max(-1.0, min(1.0, 1.0 - 2.0 * frac))))
    print(f"    Estimated critical angle: theta_crit ~ {theta_crit_est:.2f} deg")

    k_crit, sig_crit, dt_crit = run_state("bol_critical_angle", [theta_crit_est] * 10)
    print(f"    k_eff(theta_crit = {theta_crit_est:.2f} deg) = {k_crit:.5f} +/- {sig_crit:.5f} ({dt_crit:.1f} s)")

    # Export Summary CSV
    df = pd.DataFrame([{
        "T_inlet_K": 900.0,
        "T_outlet_K": 1150.0,
        "k_eff_0deg": k_0,
        "sigma_k_0deg": sig_0,
        "excess_reactivity_pcm": rho_ex_pcm,
        "sigma_excess_pcm": sig_rho_ex,
        "k_eff_180deg": k_180,
        "sigma_k_180deg": sig_180,
        "shutdown_margin_all_pcm": sdm_all_pcm,
        "sigma_sdm_all_pcm": sig_sdm_all,
        "k_eff_stuck": k_stuck,
        "sigma_k_stuck": sig_stuck,
        "shutdown_margin_stuck_pcm": sdm_stuck_pcm,
        "sigma_sdm_stuck_pcm": sig_sdm_stuck,
        "total_drum_worth_pcm": drum_worth_pcm,
        "critical_drum_angle_deg": theta_crit_est,
        "k_eff_critical": k_crit,
        "sigma_k_critical": sig_crit
    }])
    summary_path = os.path.join(DATA_DIR, "bol_neutronics_summary_third.csv")
    df.to_csv(summary_path, index=False)
    print(f"\n[+] BOL Neutronics Summary successfully exported to {summary_path}")


if __name__ == "__main__":
    main()
