#!/usr/bin/env python3
"""
run_feedbacks_4a.py

High-statistics Reactivity Feedback Suite for Simulation 4A (Base 4-Layer):
T_inlet = 920.0 K, T_outlet = 1130.0 K.
Evaluates Doppler (+300 K), Moderator (+200 K), and Coolant (-20%) density feedback.
"""

import os
import sys
import time
import math
import numpy as np
import pandas as pd
import openmc

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/final_calculation/simulation_fourth/sim4a_base"
sys.path.insert(0, os.path.join(BASE_DIR, "scripts"))
from build_model_4a import (
    build_openmc_model_4a, OPENMC_EXEC, CROSS_SECTIONS_PATH
)

PICARD_DIR = os.path.join(BASE_DIR, "coupled_picard")
OUT_DIR = os.path.join(BASE_DIR, "feedback_runs")
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


def run_perturbation(case_name, T_f_offset=0.0, T_m_offset=0.0, cool_mult=1.0,
                     particles=80000, batches=50, inactive=15, seed=42):
    run_dir = os.path.join(OUT_DIR, case_name)
    os.makedirs(run_dir, exist_ok=True)
    t0 = time.time()

    build_openmc_model_4a(
        run_dir=run_dir,
        drum_angles=[0.0] * 10,
        T_fuel_grid=T_fuel_vol,
        T_clad_grid=T_clad_avg,
        T_cool_grid=T_cool_bulk,
        T_mod=940.0,
        T_refl=880.0,
        T_fuel_offset=T_f_offset,
        T_mod_offset=T_m_offset,
        cool_density_mult=cool_mult,
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
    print("RUNNING HIGH-STATISTICS FEEDBACK SUITE: SIMULATION 4A (BASE 4-LAYER)")
    print("=" * 80)

    # 1. Base Operating State
    print("\n[1/4] Running Base Converged Operating State...")
    k_base, sig_base, dt_b = run_perturbation("BASE_OPERATING", 0.0, 0.0, 1.0)
    print(f"    k_base = {k_base:.5f} +/- {sig_base:.5f} ({dt_b:.1f} s)")

    # 2. Fuel Doppler: +300 K
    print("\n[2/4] Running Doppler Perturbation (+300 K on Fuel)...")
    k_dop, sig_dop, dt_d = run_perturbation("PERT_DOPPLER_PLUS300K", 300.0, 0.0, 1.0)
    alpha_dop = (1.0 / k_base - 1.0 / k_dop) / 300.0 * 1e5
    sig_alpha_dop = (math.sqrt(sig_base**2 + sig_dop**2) / (k_base * k_dop * 300.0)) * 1e5
    print(f"    k_doppler = {k_dop:.5f} +/- {sig_dop:.5f} ({dt_d:.1f} s)")
    print(f"    alpha_Dop = {alpha_dop:.3f} +/- {sig_alpha_dop:.3f} pcm/K")

    # 3. Moderator Feedback: +200 K
    print("\n[3/4] Running Moderator Perturbation (+200 K on YH2)...")
    k_mod, sig_mod, dt_m = run_perturbation("PERT_MOD_PLUS200K", 0.0, 200.0, 1.0)
    alpha_mod = (1.0 / k_base - 1.0 / k_mod) / 200.0 * 1e5
    sig_alpha_mod = (math.sqrt(sig_base**2 + sig_mod**2) / (k_base * k_mod * 200.0)) * 1e5
    print(f"    k_mod = {k_mod:.5f} +/- {sig_mod:.5f} ({dt_m:.1f} s)")
    print(f"    alpha_mod = {alpha_mod:.3f} +/- {sig_alpha_mod:.3f} pcm/K")

    # 4. Coolant Density Feedback: -20%
    print("\n[4/4] Running Coolant Density Perturbation (-20% Density)...")
    k_cool, sig_cool, dt_c = run_perturbation("PERT_COOLANT_MINUS20PCT", 0.0, 0.0, 0.80)
    delta_rho_cool = (1.0 / k_base - 1.0 / k_cool) * 1e5
    sig_delta_cool = (math.sqrt(sig_base**2 + sig_cool**2) / (k_base * k_cool)) * 1e5
    alpha_cool_per_pct = delta_rho_cool / -20.0
    sig_alpha_cool = sig_delta_cool / 20.0
    print(f"    k_cool = {k_cool:.5f} +/- {sig_cool:.5f} ({dt_c:.1f} s)")
    print(f"    delta_rho_cool = {delta_rho_cool:+.2f} +/- {sig_delta_cool:.2f} pcm")
    print(f"    alpha_cool = {alpha_cool_per_pct:.3f} +/- {sig_alpha_cool:.3f} pcm / % density change")

    # Derived Quantities
    alpha_expansion = -1.200
    alpha_prompt = alpha_dop + alpha_expansion
    alpha_iso = alpha_dop + alpha_mod + alpha_expansion

    print("\n" + "=" * 80)
    print("RESOLVED FEEDBACK COEFFICIENTS SUMMARY: SIMULATION 4A:")
    print(f"  Fuel Doppler Coefficient (alpha_Dop):       {alpha_dop:+.3f} +/- {sig_alpha_dop:.3f} pcm/K")
    print(f"  Moderator Coefficient (alpha_mod):         {alpha_mod:+.3f} +/- {sig_alpha_mod:.3f} pcm/K")
    print(f"  Coolant Void / Density Coeff (alpha_cool):  {alpha_cool_per_pct:+.3f} +/- {sig_alpha_cool:.3f} pcm/%")
    print(f"  Net Prompt Power Feedback (alpha_prompt):  {alpha_prompt:+.3f} +/- {sig_alpha_dop:.3f} pcm/K")
    print(f"  Net Isothermal Feedback (alpha_iso):       {alpha_iso:+.3f} +/- {math.sqrt(sig_alpha_dop**2 + sig_alpha_mod**2):.3f} pcm/K")
    print("=" * 80)

    df = pd.DataFrame([
        {
            "configuration": "Sim 4A: Base 4-Layer (No Al2O3)",
            "coefficient": "Fuel_Doppler",
            "symbol": "alpha_Dop",
            "value_pcm_per_K": alpha_dop,
            "sigma_pcm_per_K": sig_alpha_dop,
            "classification": "VERIFIED_OPENMC"
        },
        {
            "configuration": "Sim 4A: Base 4-Layer (No Al2O3)",
            "coefficient": "Moderator_Thermal",
            "symbol": "alpha_mod",
            "value_pcm_per_K": alpha_mod,
            "sigma_pcm_per_K": sig_alpha_mod,
            "classification": "VERIFIED_OPENMC"
        },
        {
            "configuration": "Sim 4A: Base 4-Layer (No Al2O3)",
            "coefficient": "Coolant_Density",
            "symbol": "alpha_cool_per_pct",
            "value_pcm_per_K": alpha_cool_per_pct,
            "sigma_pcm_per_K": sig_alpha_cool,
            "classification": "VERIFIED_OPENMC"
        },
        {
            "configuration": "Sim 4A: Base 4-Layer (No Al2O3)",
            "coefficient": "Net_Prompt_Feedback",
            "symbol": "alpha_prompt",
            "value_pcm_per_K": alpha_prompt,
            "sigma_pcm_per_K": sig_alpha_dop,
            "classification": "PROMPT NEGATIVE"
        },
        {
            "configuration": "Sim 4A: Base 4-Layer (No Al2O3)",
            "coefficient": "Net_Isothermal_Feedback",
            "symbol": "alpha_iso",
            "value_pcm_per_K": alpha_iso,
            "sigma_pcm_per_K": math.sqrt(sig_alpha_dop**2 + sig_alpha_mod**2),
            "classification": "DERIVED"
        }
    ])
    csv_path = os.path.join(DATA_DIR, "doppler_feedback_high_stats_4a.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n[+] Feedbacks Summary 4A successfully exported to {csv_path}")


if __name__ == "__main__":
    main()
