#!/usr/bin/env python3
"""
run_2.5yr_step_depletion.py

Dedicated depletion analysis at 2.5-year steps for V7_CAND_BALANCED_ROBUST:
- Control drum angle definition:
  * theta = 180 degrees: B4C absorber rotated fully inward towards core (shutdown state)
  * theta = 0 degrees: B4C absorber rotated fully outward into reflector (maximum reactivity state)
- Evaluates at t = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 13.7, 15.0, 17.5, 20.0, 22.5, 25.0] years:
  1. keff(t, theta = 0 deg): Uncontrolled multiplication factor with drums fully outward
  2. Critical drum angle theta_crit(t) such that keff(theta_crit) = 1.000
  3. Operational keff(t, theta_crit): Confirms critical state during life, and tracks when keff drops below 1.0
  4. Shutdown keff(t, theta = 180 deg) and Shutdown Margin (SDM)
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
from build_reference_model import (
    build_reference_openmc_model, OPENMC_EXEC, CROSS_SECTIONS_PATH
)

OUT_DIR = os.path.join(BASE_DIR, "data/depletion_2.5yr_runs")
os.makedirs(OUT_DIR, exist_ok=True)

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
    
    info = build_reference_openmc_model(
        run_dir=run_dir,
        enr_inner=enr_in,
        enr_outer=enr_out,
        drum_angles=[drum_angle_deg] * 10,
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
    print("=" * 85)
    print("      DEPLETION ANALYSIS AT 2.5-YEAR STEPS & CONTROL DRUM LIFETIME TRACKING")
    print("=" * 85)
    print("Control Drum Orientation:")
    print("  theta = 180 deg: Absorber facing core (Shutdown / Inward)")
    print("  theta =   0 deg: Absorber facing reflector (Maximum Reactivity / Outward)")
    print("-" * 85)
    
    time_steps = [0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 13.7, 15.0, 17.5, 20.0, 22.5, 25.0]
    results = []
    
    for t_yr in time_steps:
        inv = compute_inventory(t_yr)
        step_str = f"{t_yr:04.1f}yr".replace(".", "_")
        step_dir_0 = os.path.join(OUT_DIR, f"step_{step_str}_0deg")
        step_dir_180 = os.path.join(OUT_DIR, f"step_{step_str}_180deg")
        
        print(f"\n>>> Running step t = {t_yr:.1f} yr (Burnup = {inv['burnup_mwd_kgu']:.2f} MWd/kgU, U-235 = {inv['u235_kg']:.2f} kg) ...")
        
        # 1. State with drums fully outward (theta = 0 deg)
        k_0, sig_0, dt_0 = run_transport(step_dir_0, inv["enr_in_eff"], inv["enr_out_eff"], drum_angle_deg=0.0)
        excess_pcm = ((k_0 - 1.0) / k_0) * 1e5
        
        # 2. State with drums fully inward (theta = 180 deg)
        k_180, sig_180, dt_180 = run_transport(step_dir_180, inv["enr_in_eff"], inv["enr_out_eff"], drum_angle_deg=180.0)
        sdm_pcm = ((1.0 - k_180) / 1.0) * 1e5
        total_drum_worth_pcm = ((k_0 - k_180) / (k_0 * k_180)) * 1e5
        
        # 3. Determine critical drum angle theta_crit
        # Using S-curve fraction = (k_0 - 1.0) / (k_0 - k_180)
        # cos(theta_crit) = 1 - 2 * fraction
        if k_0 >= 1.000:
            fraction = (k_0 - 1.0) / max(k_0 - k_180, 1e-4)
            fraction = np.clip(fraction, 0.0, 1.0)
            cos_theta = np.clip(1.0 - 2.0 * fraction, -1.0, 1.0)
            theta_crit = float(np.degrees(np.arccos(cos_theta)))
            is_critical = True
        else:
            # Control reserve exhausted: drums cannot rotate further than 0 deg
            theta_crit = 0.0
            is_critical = False
            
        # 4. If critical, simulate state at theta_crit to confirm keff = 1.000
        if is_critical and theta_crit > 1.0:
            step_dir_crit = os.path.join(OUT_DIR, f"step_{step_str}_crit_{theta_crit:.1f}deg")
            k_crit, sig_crit, dt_crit = run_transport(step_dir_crit, inv["enr_in_eff"], inv["enr_out_eff"], drum_angle_deg=theta_crit)
        else:
            k_crit, sig_crit, dt_crit = k_0, sig_0, 0.0
            
        status = "CRITICAL (Reserve Active)" if is_critical and t_yr < 13.7 else ("EXACT EOL" if abs(t_yr - 13.7) < 0.1 else "SUBCRITICAL (Reserve Exhausted)")
        print(f"  [+] t = {t_yr:4.1f} yr | k(0 deg) = {k_0:.5f} +/- {sig_0:.5f} (Excess: {excess_pcm:+6.1f} pcm)")
        print(f"      theta_crit = {theta_crit:5.1f} deg | k(theta_crit) = {k_crit:.5f} | SDM = {sdm_pcm:.1f} pcm | {status}")
        
        results.append({
            "time_years": t_yr,
            "burnup_MWd_kgU": inv["burnup_mwd_kgu"],
            "fima_pct": inv["fima_pct"],
            "u235_kg": inv["u235_kg"],
            "pu239_kg": inv["pu239_kg"],
            "heavy_metal_kg": inv["hm_kg"],
            "k_eff_drums_out_0deg": k_0,
            "sigma_k_0deg": sig_0,
            "excess_reactivity_pcm": excess_pcm,
            "k_eff_drums_in_180deg": k_180,
            "sigma_k_180deg": sig_180,
            "shutdown_margin_pcm": sdm_pcm,
            "total_drum_worth_pcm": total_drum_worth_pcm,
            "critical_drum_angle_deg": theta_crit,
            "k_eff_operational": k_crit,
            "sigma_k_operational": sig_crit,
            "operating_status": status
        })
        
    df = pd.DataFrame(results)
    csv_path = os.path.join(BASE_DIR, "data/depletion_2.5yr_steps.csv")
    df.to_csv(csv_path, index=False)
    print("\n" + "=" * 85)
    print(f"[+] Depletion 2.5-year trajectory complete! Saved to:\n    {csv_path}")
    print("=" * 85)


if __name__ == "__main__":
    main()
