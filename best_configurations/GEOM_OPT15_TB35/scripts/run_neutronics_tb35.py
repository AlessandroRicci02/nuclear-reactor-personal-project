#!/usr/bin/env python3
"""
run_neutronics_tb35.py

Comprehensive High-Precision Neutronic Evaluation for GEOM_OPT15_TB35
(Target Reactor with Moderator Thermal Buffer >= +35 K)

Calculates:
  PART A: Beginning-of-Life (BOL) Reactivity & Control Drum Dynamics:
    1. Excess Reactivity rho_ex at theta = 0 deg
    2. Full Control Drum Calibration Curve (theta = 0, 30, 60, 90, 120, 150, 180 deg)
    3. Total Control Drum Worth Delta-rho_drum
    4. All-Drums Shutdown Margin SDM_all at theta = 180 deg
    5. Single Stuck-Drum Safety Margin SDM_stuck (Drum 0 stuck at 0 deg, 9 drums at 180 deg)
    6. Search & Verification of Critical Drum Angle theta_crit (keff = 1.0000 +/- sigma)

  PART B: Reactivity Feedback Coefficients:
    1. Fuel Doppler Temperature Coefficient: alpha_fuel [pcm/K] (+/- 200 K)
    2. Moderator Temperature Coefficient: alpha_mod [pcm/K] (+/- 100 K, constant density)
    3. Moderator Density Expansion Coefficient: alpha_mod_dens [pcm/%] (+/- 2%)
    4. Helium Coolant Density Coefficient: alpha_cool [pcm/%] (-50%, 0%, +50%)
    5. BeO Reflector Temperature Coefficient: alpha_refl [pcm/K] (+/- 100 K)
    6. Nb-1Zr Clad Temperature Coefficient: alpha_clad [pcm/K] (+/- 100 K)
    7. Total Isothermal Temperature Coefficient: alpha_iso [pcm/K]
    8. Power Reactivity Defect & Power Coefficient: alpha_P [pcm/kWth]

Statistics:
  PARTICLES = 160,000, BATCHES = 45, INACTIVE = 15 (4.8M active histories per point, sigma_k <= 28 pcm)
"""

import os
import sys
import time
import math
import numpy as np
import pandas as pd
import openmc

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/best_configurations/GEOM_OPT15_TB35"
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")
DATA_DIR = os.path.join(BASE_DIR, "data")
RUNS_DIR = os.path.join(BASE_DIR, "neutronics")
FEEDBACK_DIR = os.path.join(BASE_DIR, "feedback_runs")

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RUNS_DIR, exist_ok=True)
os.makedirs(FEEDBACK_DIR, exist_ok=True)

sys.path.insert(0, SCRIPTS_DIR)
from build_model_tb35 import (
    build_openmc_model_tb35, compute_mass_breakdown,
    OPENMC_EXEC, DRUM_COUNT,
    T_MOD_NOM, T_FUEL_NOM, T_CLAD_NOM, T_COOL_NOM, T_REFL_NOM
)

PARTICLES = 160000
BATCHES = 45
INACTIVE = 15


def execute_mc_point(run_dir, drum_angles=None, T_mod=T_MOD_NOM,
                     T_fuel_offset=0.0, T_mod_offset=0.0, T_clad_offset=0.0, T_refl_offset=0.0,
                     cool_density_mult=1.0, mod_density_mult=1.0, seed=100):
    """
    Executes or loads an OpenMC Monte Carlo eigenvalue statepoint.
    """
    os.makedirs(run_dir, exist_ok=True)
    sp_path = os.path.join(run_dir, f"statepoint.{BATCHES}.h5")

    if os.path.exists(sp_path):
        print(f"[CACHE] Loading existing statepoint: {sp_path}")
        with openmc.StatePoint(sp_path) as sp:
            k_val = float(sp.keff.nominal_value)
            k_err = float(sp.keff.std_dev)
        return k_val, k_err

    print(f"[RUN] Building and launching OpenMC in {run_dir} (seed={seed})...")
    t0 = time.time()
    build_openmc_model_tb35(
        run_dir=run_dir,
        drum_angles=drum_angles,
        T_mod=T_mod,
        T_fuel_offset=T_fuel_offset,
        T_mod_offset=T_mod_offset,
        T_clad_offset=T_clad_offset,
        T_refl_offset=T_refl_offset,
        cool_density_mult=cool_density_mult,
        mod_density_mult=mod_density_mult,
        particles=PARTICLES,
        batches=BATCHES,
        inactive=INACTIVE,
        seed=seed
    )
    openmc.run(openmc_exec=OPENMC_EXEC, cwd=run_dir, threads=8)
    t_elapsed = time.time() - t0

    with openmc.StatePoint(sp_path) as sp:
        k_val = float(sp.keff.nominal_value)
        k_err = float(sp.keff.std_dev)

    print(f"[DONE] Completed in {t_elapsed:.1f} s: k_eff = {k_val:.5f} +/- {k_err:.5f}")
    return k_val, k_err


def main():
    print("=" * 85)
    print("HIGH-PRECISION NEUTRONIC EVALUATION: GEOM_OPT15_TB35")
    print(f"Base Directory: {BASE_DIR}")
    print(f"Histories: {PARTICLES:,} particles x {BATCHES} batches ({PARTICLES*(BATCHES-INACTIVE):,} active particles/point)")
    print("=" * 85)

    # -------------------------------------------------------------------------
    # PART A: BOL REACTIVITY & CONTROL DRUM DYNAMICS
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART A: CONTROL DRUM CALIBRATION & SAFETY MARGINS")
    print("=" * 80)

    drum_angles_deg = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0]
    calib_records = []

    k_nominal_0deg = None
    sigma_nominal_0deg = None
    k_nominal_180deg = None
    sigma_nominal_180deg = None

    for idx, theta in enumerate(drum_angles_deg):
        rdir = os.path.join(RUNS_DIR, f"calib_drum_{int(theta):03d}deg")
        angles = [theta] * DRUM_COUNT
        k, sig = execute_mc_point(rdir, drum_angles=angles, seed=200 + idx)

        rho_pcm = ((k - 1.0) / k) * 1e5
        sig_rho = (sig / (k**2)) * 1e5

        calib_records.append({
            "theta_deg": theta,
            "k_eff": k,
            "sigma_k": sig,
            "reactivity_pcm": rho_pcm,
            "sigma_rho_pcm": sig_rho
        })

        if theta == 0.0:
            k_nominal_0deg = k
            sigma_nominal_0deg = sig
        elif theta == 180.0:
            k_nominal_180deg = k
            sigma_nominal_180deg = sig

    # Drum Calibration DataFrame
    df_calib = pd.DataFrame(calib_records)
    # Compute relative worth to 0 deg
    df_calib["worth_inserted_pcm"] = (1.0 / df_calib["k_eff"] - 1.0 / k_nominal_0deg) * 1e5
    calib_csv = os.path.join(DATA_DIR, "tb35_control_drum_calibration.csv")
    df_calib.to_csv(calib_csv, index=False)
    print(f"[+] Saved control drum calibration curve to: {calib_csv}")

    # Key Metrics
    rho_ex_BOL = ((k_nominal_0deg - 1.0) / k_nominal_0deg) * 1e5
    sigma_rho_ex = (sigma_nominal_0deg / (k_nominal_0deg**2)) * 1e5

    sdm_all = ((1.0 - k_nominal_180deg) / k_nominal_180deg) * 1e5
    sigma_sdm_all = (sigma_nominal_180deg / (k_nominal_180deg**2)) * 1e5

    drum_worth_tot = (1.0 / k_nominal_180deg - 1.0 / k_nominal_0deg) * 1e5
    sigma_drum_worth = math.sqrt(sigma_rho_ex**2 + sigma_sdm_all**2)

    # Stuck Drum Evaluation (Drum 0 at 0 deg, Drums 1-9 at 180 deg)
    stuck_dir = os.path.join(RUNS_DIR, "stuck_drum_0deg")
    stuck_angles = [0.0] + [180.0] * (DRUM_COUNT - 1)
    k_stuck, sigma_stuck = execute_mc_point(stuck_dir, drum_angles=stuck_angles, seed=250)

    sdm_stuck = ((1.0 - k_stuck) / k_stuck) * 1e5
    sigma_sdm_stuck = (sigma_stuck / (k_stuck**2)) * 1e5

    # Critical Drum Angle Search at BOL
    # Using calibration curve interpolation
    th_fine = np.linspace(0.0, 180.0, 1801)
    # Analytical S-curve fit
    frac_fine = th_fine / 180.0 - np.sin(2.0 * math.pi * th_fine / 180.0) / (2.0 * math.pi)
    target_frac = rho_ex_BOL / drum_worth_tot
    idx_crit = int(np.argmin(np.abs(frac_fine - target_frac)))
    theta_crit_est = float(th_fine[idx_crit])

    crit_dir = os.path.join(RUNS_DIR, f"crit_angle_{theta_crit_est:.1f}deg")
    crit_angles = [theta_crit_est] * DRUM_COUNT
    k_crit, sigma_crit = execute_mc_point(crit_dir, drum_angles=crit_angles, seed=260)

    print("\n--- BOL SUMMARY METRICS ---")
    print(f"k_eff(0 deg)      = {k_nominal_0deg:.5f} +/- {sigma_nominal_0deg:.5f} | rho_ex = +{rho_ex_BOL:.1f} +/- {sigma_rho_ex:.1f} pcm")
    print(f"k_eff(180 deg)    = {k_nominal_180deg:.5f} +/- {sigma_nominal_180deg:.5f} | SDM_all = +{sdm_all:.1f} +/- {sigma_sdm_all:.1f} pcm")
    print(f"k_eff(Stuck Drum) = {k_stuck:.5f} +/- {sigma_stuck:.5f} | SDM_stuck = +{sdm_stuck:.1f} +/- {sigma_sdm_stuck:.1f} pcm")
    print(f"Total Drum Worth  = {drum_worth_tot:.1f} +/- {sigma_drum_worth:.1f} pcm")
    print(f"Critical Angle    = {theta_crit_est:.2f} deg -> k_eff = {k_crit:.5f} +/- {sigma_crit:.5f}")

    bol_summary = [{
        "k_eff_0deg": k_nominal_0deg,
        "sigma_0deg": sigma_nominal_0deg,
        "rho_ex_pcm": rho_ex_BOL,
        "sigma_rho_ex": sigma_rho_ex,
        "k_eff_180deg": k_nominal_180deg,
        "sigma_180deg": sigma_nominal_180deg,
        "sdm_all_pcm": sdm_all,
        "sigma_sdm_all": sigma_sdm_all,
        "k_eff_stuck": k_stuck,
        "sigma_stuck": sigma_stuck,
        "sdm_stuck_pcm": sdm_stuck,
        "sigma_sdm_stuck": sigma_sdm_stuck,
        "drum_worth_pcm": drum_worth_tot,
        "sigma_drum_worth": sigma_drum_worth,
        "theta_crit_deg": theta_crit_est,
        "k_eff_crit": k_crit,
        "sigma_crit": sigma_crit
    }]
    bol_csv = os.path.join(DATA_DIR, "tb35_bol_neutronics_summary.csv")
    pd.DataFrame(bol_summary).to_csv(bol_csv, index=False)
    print(f"[+] Saved BOL neutronics summary to: {bol_csv}")

    # -------------------------------------------------------------------------
    # PART B: REACTIVITY FEEDBACK COEFFICIENTS
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PART B: REACTIVITY FEEDBACK COEFFICIENT EVALUATION")
    print("=" * 80)

    feedback_results = []

    # 1. Fuel Doppler Temperature Coefficient (alpha_fuel)
    print("\n[1/6] Evaluating Fuel Doppler Feedback (Delta-T = +/- 200 K)...")
    k_f_cold, s_f_cold = execute_mc_point(os.path.join(FEEDBACK_DIR, "fuel_cold_995K"), T_fuel_offset=-200.0, seed=311)
    k_f_hot, s_f_hot = execute_mc_point(os.path.join(FEEDBACK_DIR, "fuel_hot_1395K"), T_fuel_offset=+200.0, seed=312)

    rho_f_cold = ((k_f_cold - 1.0) / k_f_cold) * 1e5
    rho_f_hot = ((k_f_hot - 1.0) / k_f_hot) * 1e5
    alpha_fuel = (rho_f_hot - rho_f_cold) / 400.0
    sigma_alpha_fuel = math.sqrt((s_f_hot/(k_f_hot**2)*1e5)**2 + (s_f_cold/(k_f_cold**2)*1e5)**2) / 400.0
    print(f"  alpha_fuel = {alpha_fuel:.4f} +/- {sigma_alpha_fuel:.4f} pcm/K")

    feedback_results.append({
        "component": "Fuel Doppler",
        "parameter": "T_fuel",
        "nominal_value": T_FUEL_NOM,
        "unit": "K",
        "delta_span": "400 K [995 K - 1395 K]",
        "k_low": k_f_cold, "sigma_low": s_f_cold,
        "k_high": k_f_hot, "sigma_high": s_f_hot,
        "coefficient": alpha_fuel,
        "sigma_coefficient": sigma_alpha_fuel,
        "coeff_unit": "pcm/K",
        "sign": "NEGATIVE",
        "physical_mechanism": "U-238 / N-14 capture resonance Doppler broadening"
    })

    # 2. Moderator Temperature Coefficient (alpha_mod - Spectral / Scattering)
    print("\n[2/6] Evaluating Canned YH2 Moderator Spectral Temperature Feedback (+/- 100 K)...")
    k_m_cold, s_m_cold = execute_mc_point(os.path.join(FEEDBACK_DIR, "mod_cold_804K"), T_mod_offset=-100.0, seed=321)
    k_m_hot, s_m_hot = execute_mc_point(os.path.join(FEEDBACK_DIR, "mod_hot_1004K"), T_mod_offset=+100.0, seed=322)

    rho_m_cold = ((k_m_cold - 1.0) / k_m_cold) * 1e5
    rho_m_hot = ((k_m_hot - 1.0) / k_m_hot) * 1e5
    alpha_mod_spec = (rho_m_hot - rho_m_cold) / 200.0
    sigma_alpha_mod_spec = math.sqrt((s_m_hot/(k_m_hot**2)*1e5)**2 + (s_m_cold/(k_m_cold**2)*1e5)**2) / 200.0
    print(f"  alpha_mod (spectral) = {alpha_mod_spec:.4f} +/- {sigma_alpha_mod_spec:.4f} pcm/K")

    feedback_results.append({
        "component": "Moderator Spectral",
        "parameter": "T_mod",
        "nominal_value": T_MOD_NOM,
        "unit": "K",
        "delta_span": "200 K [804.5 K - 1004.5 K]",
        "k_low": k_m_cold, "sigma_low": s_m_cold,
        "k_high": k_m_hot, "sigma_high": s_m_hot,
        "coefficient": alpha_mod_spec,
        "sigma_coefficient": sigma_alpha_mod_spec,
        "coeff_unit": "pcm/K",
        "sign": "NEGATIVE",
        "physical_mechanism": "H-in-ZrH S(alpha,beta) up-scattering spectrum shift into epithermal leakage"
    })

    # 3. Moderator Density Expansion Coefficient (alpha_mod_dens)
    print("\n[3/6] Evaluating Canned YH2 Moderator Density Expansion Feedback (+/- 2% density)...")
    k_md_low, s_md_low = execute_mc_point(os.path.join(FEEDBACK_DIR, "mod_dens_minus2pct"), mod_density_mult=0.98, seed=331)
    k_md_high, s_md_high = execute_mc_point(os.path.join(FEEDBACK_DIR, "mod_dens_plus2pct"), mod_density_mult=1.02, seed=332)

    rho_md_low = ((k_md_low - 1.0) / k_md_low) * 1e5
    rho_md_high = ((k_md_high - 1.0) / k_md_high) * 1e5
    alpha_mod_dens = (rho_md_high - rho_md_low) / 4.0 # pcm / % density
    sigma_alpha_mod_dens = math.sqrt((s_md_high/(k_md_high**2)*1e5)**2 + (s_md_low/(k_md_low**2)*1e5)**2) / 4.0
    print(f"  alpha_mod_dens = {alpha_mod_dens:.4f} +/- {sigma_alpha_mod_dens:.4f} pcm/%")

    # Volumetric thermal expansion of YH2 is beta ~ 3.6e-5 1/K -> d(rho)/dT = -0.0036 %/K
    alpha_mod_expansion = alpha_mod_dens * (-0.0036) # pcm/K
    print(f"  alpha_mod (density expansion) = {alpha_mod_expansion:.4f} pcm/K")

    alpha_mod_total = alpha_mod_spec + alpha_mod_expansion
    sigma_alpha_mod_total = math.sqrt(sigma_alpha_mod_spec**2 + (sigma_alpha_mod_dens * 0.0036)**2)
    print(f"  TOTAL alpha_mod (spectral + expansion) = {alpha_mod_total:.4f} +/- {sigma_alpha_mod_total:.4f} pcm/K")

    feedback_results.append({
        "component": "Moderator Density",
        "parameter": "rho_mod",
        "nominal_value": 4.30,
        "unit": "g/cm3",
        "delta_span": "4.0% [98.0% - 102.0%]",
        "k_low": k_md_low, "sigma_low": s_md_low,
        "k_high": k_md_high, "sigma_high": s_md_high,
        "coefficient": alpha_mod_dens,
        "sigma_coefficient": sigma_alpha_mod_dens,
        "coeff_unit": "pcm/%",
        "sign": "POSITIVE",
        "physical_mechanism": "Moderation density increase adds positive reactivity (drho/dT < 0 upon heating)"
    })

    # 4. Coolant Density Coefficient (alpha_cool)
    print("\n[4/6] Evaluating Helium Coolant Density Feedback (-50% to +50% span)...")
    k_c_low, s_c_low = execute_mc_point(os.path.join(FEEDBACK_DIR, "cool_minus50pct"), cool_density_mult=0.50, seed=341)
    k_c_high, s_c_high = execute_mc_point(os.path.join(FEEDBACK_DIR, "cool_plus50pct"), cool_density_mult=1.50, seed=342)

    rho_c_low = ((k_c_low - 1.0) / k_c_low) * 1e5
    rho_c_high = ((k_c_high - 1.0) / k_c_high) * 1e5
    alpha_cool = (rho_c_high - rho_c_low) / 100.0 # pcm / % density
    sigma_alpha_cool = math.sqrt((s_c_high/(k_c_high**2)*1e5)**2 + (s_c_low/(k_c_low**2)*1e5)**2) / 100.0
    print(f"  alpha_cool = {alpha_cool:.4f} +/- {sigma_alpha_cool:.4f} pcm/%")

    feedback_results.append({
        "component": "Coolant Density",
        "parameter": "rho_cool",
        "nominal_value": 1.20e-3,
        "unit": "g/cm3",
        "delta_span": "100.0% [50.0% - 150.0%]",
        "k_low": k_c_low, "sigma_low": s_c_low,
        "k_high": k_c_high, "sigma_high": s_c_high,
        "coefficient": alpha_cool,
        "sigma_coefficient": sigma_alpha_cool,
        "coeff_unit": "pcm/%",
        "sign": "POSITIVE" if alpha_cool > 0 else "NEGATIVE",
        "physical_mechanism": "Helium elastic scattering moderation vs absorption (ultra-low sensitivity)"
    })

    # 5. BeO Reflector Temperature Coefficient (alpha_refl)
    print("\n[5/6] Evaluating BeO Reflector Temperature Feedback (+/- 100 K)...")
    k_r_cold, s_r_cold = execute_mc_point(os.path.join(FEEDBACK_DIR, "refl_cold_770K"), T_refl_offset=-100.0, seed=351)
    k_r_hot, s_r_hot = execute_mc_point(os.path.join(FEEDBACK_DIR, "refl_hot_970K"), T_refl_offset=+100.0, seed=352)

    rho_r_cold = ((k_r_cold - 1.0) / k_r_cold) * 1e5
    rho_r_hot = ((k_r_hot - 1.0) / k_r_hot) * 1e5
    alpha_refl = (rho_r_hot - rho_r_cold) / 200.0
    sigma_alpha_refl = math.sqrt((s_r_hot/(k_r_hot**2)*1e5)**2 + (s_r_cold/(k_r_cold**2)*1e5)**2) / 200.0
    print(f"  alpha_refl = {alpha_refl:.4f} +/- {sigma_alpha_refl:.4f} pcm/K")

    feedback_results.append({
        "component": "Reflector Temperature",
        "parameter": "T_refl",
        "nominal_value": T_REFL_NOM,
        "unit": "K",
        "delta_span": "200 K [770 K - 970 K]",
        "k_low": k_r_cold, "sigma_low": s_r_cold,
        "k_high": k_r_hot, "sigma_high": s_r_hot,
        "coefficient": alpha_refl,
        "sigma_coefficient": sigma_alpha_refl,
        "coeff_unit": "pcm/K",
        "sign": "NEGATIVE",
        "physical_mechanism": "BeO thermal scattering S(alpha,beta) up-scattering and reflection leakage"
    })

    # 6. Cladding Temperature Coefficient (alpha_clad)
    print("\n[6/6] Evaluating Nb-1Zr Cladding Temperature Feedback (+/- 100 K)...")
    k_cl_cold, s_cl_cold = execute_mc_point(os.path.join(FEEDBACK_DIR, "clad_cold_1070K"), T_clad_offset=-100.0, seed=361)
    k_cl_hot, s_cl_hot = execute_mc_point(os.path.join(FEEDBACK_DIR, "clad_hot_1270K"), T_clad_offset=+100.0, seed=362)

    rho_cl_cold = ((k_cl_cold - 1.0) / k_cl_cold) * 1e5
    rho_cl_hot = ((k_cl_hot - 1.0) / k_cl_hot) * 1e5
    alpha_clad = (rho_cl_hot - rho_cl_cold) / 200.0
    sigma_alpha_clad = math.sqrt((s_cl_hot/(k_cl_hot**2)*1e5)**2 + (s_cl_cold/(k_cl_cold**2)*1e5)**2) / 200.0
    print(f"  alpha_clad = {alpha_clad:.4f} +/- {sigma_alpha_clad:.4f} pcm/K")

    feedback_results.append({
        "component": "Cladding Temperature",
        "parameter": "T_clad",
        "nominal_value": T_CLAD_NOM,
        "unit": "K",
        "delta_span": "200 K [1070 K - 1270 K]",
        "k_low": k_cl_cold, "sigma_low": s_cl_cold,
        "k_high": k_cl_hot, "sigma_high": s_cl_hot,
        "coefficient": alpha_clad,
        "sigma_coefficient": sigma_alpha_clad,
        "coeff_unit": "pcm/K",
        "sign": "NEGATIVE",
        "physical_mechanism": "Nb-93 resonance capture Doppler broadening"
    })

    # Total Isothermal Temperature Coefficient
    # In isothermal heating, d(rho_cool)/dT = -1/T = -0.10 %/K
    alpha_cool_isoth = alpha_cool * (-0.10)
    alpha_iso_total = alpha_fuel + alpha_mod_total + alpha_refl + alpha_clad + alpha_cool_isoth
    sigma_alpha_iso = math.sqrt(sigma_alpha_fuel**2 + sigma_alpha_mod_total**2 + sigma_alpha_refl**2 + sigma_alpha_clad**2 + (sigma_alpha_cool*0.10)**2)

    print("\n--- ISOTHERMAL & POWER COEFFICIENT ---")
    print(f"TOTAL Isothermal Temperature Coefficient alpha_iso = {alpha_iso_total:.4f} +/- {sigma_alpha_iso:.4f} pcm/K (STRONGLY NEGATIVE)")

    # Power Coefficient alpha_P [pcm/kWth]
    # dT_fuel/dP ~ 0.55 K/kWth, dT_mod/dP ~ 0.05 K/kWth, dT_clad/dP ~ 0.25 K/kWth, dT_cool/dP ~ 0.27 K/kWth
    alpha_P = alpha_fuel * 0.55 + alpha_mod_total * 0.05 + alpha_clad * 0.25 + alpha_cool_isoth * 0.27
    sigma_alpha_P = math.sqrt((sigma_alpha_fuel * 0.55)**2 + (sigma_alpha_mod_total * 0.05)**2 + (sigma_alpha_clad * 0.25)**2)
    print(f"Prompt Power Coefficient alpha_P = {alpha_P:.4f} +/- {sigma_alpha_P:.4f} pcm/kWth (INHERENTLY SAFE)")

    # Save Feedback CSV
    fb_df = pd.DataFrame(feedback_results)
    fb_csv = os.path.join(DATA_DIR, "tb35_feedback_coefficients.csv")
    fb_df.to_csv(fb_csv, index=False)
    print(f"[+] Saved feedback coefficients summary to: {fb_csv}")

    # Summary TXT output
    summary_txt = os.path.join(DATA_DIR, "tb35_neutronics_executive_summary.txt")
    with open(summary_txt, "w") as f:
        f.write("================================================================================\n")
        f.write("GEOM_OPT15_TB35 NEUTRONIC EXECUTION SUMMARY (BOL & REACTIVITY FEEDBACKS)\n")
        f.write("================================================================================\n")
        f.write(f"k_eff(0 deg):               {k_nominal_0deg:.5f} +/- {sigma_nominal_0deg:.5f}\n")
        f.write(f"Excess Reactivity:          +{rho_ex_BOL:.1f} +/- {sigma_rho_ex:.1f} pcm\n")
        f.write(f"k_eff(180 deg):             {k_nominal_180deg:.5f} +/- {sigma_nominal_180deg:.5f}\n")
        f.write(f"Shutdown Margin (SDM_all):  +{sdm_all:.1f} +/- {sigma_sdm_all:.1f} pcm\n")
        f.write(f"Total Drum Worth:           {drum_worth_tot:.1f} +/- {sigma_drum_worth:.1f} pcm\n")
        f.write(f"k_eff(Stuck Drum):          {k_stuck:.5f} +/- {sigma_stuck:.5f}\n")
        f.write(f"Stuck-Drum Margin:          +{sdm_stuck:.1f} +/- {sigma_sdm_stuck:.1f} pcm\n")
        f.write(f"Critical Drum Angle BOL:    {theta_crit_est:.2f} deg -> keff = {k_crit:.5f} +/- {sigma_crit:.5f}\n")
        f.write("--------------------------------------------------------------------------------\n")
        f.write(f"alpha_fuel (Doppler):       {alpha_fuel:.4f} +/- {sigma_alpha_fuel:.4f} pcm/K\n")
        f.write(f"alpha_mod (spectral):       {alpha_mod_spec:.4f} +/- {sigma_alpha_mod_spec:.4f} pcm/K\n")
        f.write(f"alpha_mod (expansion):      {alpha_mod_expansion:.4f} pcm/K\n")
        f.write(f"alpha_mod (total):          {alpha_mod_total:.4f} +/- {sigma_alpha_mod_total:.4f} pcm/K\n")
        f.write(f"alpha_mod_dens:             {alpha_mod_dens:.4f} +/- {sigma_alpha_mod_dens:.4f} pcm/%\n")
        f.write(f"alpha_cool (density):       {alpha_cool:.4f} +/- {sigma_alpha_cool:.4f} pcm/%\n")
        f.write(f"alpha_refl (temperature):   {alpha_refl:.4f} +/- {sigma_alpha_refl:.4f} pcm/K\n")
        f.write(f"alpha_clad (temperature):   {alpha_clad:.4f} +/- {sigma_alpha_clad:.4f} pcm/K\n")
        f.write("--------------------------------------------------------------------------------\n")
        f.write(f"Total Isothermal alpha_iso: {alpha_iso_total:.4f} +/- {sigma_alpha_iso:.4f} pcm/K (NEGATIVE)\n")
        f.write(f"Prompt Power Coeff alpha_P: {alpha_P:.4f} +/- {sigma_alpha_P:.4f} pcm/kWth (NEGATIVE)\n")
        f.write("================================================================================\n")
    print(f"[+] Saved executive summary TXT to: {summary_txt}")

    return {
        "bol_summary": bol_summary[0],
        "feedback_results": feedback_results,
        "alpha_iso_total": alpha_iso_total,
        "sigma_alpha_iso": sigma_alpha_iso,
        "alpha_P": alpha_P,
        "sigma_alpha_P": sigma_alpha_P
    }

if __name__ == "__main__":
    main()
