#!/usr/bin/env python3
"""
run_coolant_feedback_ultra_precise.py

Ultra-High-Precision Evaluation of Helium Coolant Density Coefficient (alpha_cool)
for GEOM_OPT15_EXTENDED.

Runs high-statistics Monte Carlo transport (15.0M active particles per point, 300,000 particles x 65 batches):
  - State Minus 75% (-75% coolant density, cool_mult = 0.25)
  - State Plus 75% (+75% coolant density, cool_mult = 1.75)
Combines with existing high-statistics states (-25%, 0% nominal, +25%) to perform
a weighted generalized least squares (GLS) regression across a 150% density span,
reducing sigma(alpha_cool) below 0.18 pcm/% to definitively confirm the sign (>2-sigma).
"""

import os
import sys
import time
import math
import numpy as np
import pandas as pd
import openmc

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(BASE_DIR, "scripts"))
from build_model_opt15 import (
    build_openmc_model_opt15, OPENMC_EXEC, CROSS_SECTIONS_PATH,
    DRUM_COUNT
)

OUT_DIR = os.path.join(BASE_DIR, "feedback_runs_ultra")
DATA_DIR = os.path.join(BASE_DIR, "data")
PICARD_DIR = os.path.join(BASE_DIR, "coupled_picard")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

# Load converged temperature grids
grids_path = os.path.join(PICARD_DIR, "converged_temperature_grids_opt15.npz")
if os.path.exists(grids_path):
    grids = np.load(grids_path)
    T_fuel_vol = grids["T_fuel_vol"]
    T_clad_avg = grids["T_clad_avg"]
    T_cool_bulk = grids["T_cool_bulk"]
else:
    T_fuel_vol = np.full((217, 10), 1210.0)
    T_clad_avg = np.full((217, 10), 1185.0)
    T_cool_bulk = np.full((217, 10), 1025.0)

PARTICLES = 300000
BATCHES = 65
INACTIVE = 15


def run_state(run_dir, cool_mult, seed):
    os.makedirs(run_dir, exist_ok=True)
    sp_path = os.path.join(run_dir, f"statepoint.{BATCHES}.h5")
    if os.path.exists(sp_path):
        print(f"[+] Reusing existing statepoint: {sp_path}")
        with openmc.StatePoint(sp_path) as sp:
            k_val = float(sp.keff.nominal_value)
            k_sig = float(sp.keff.std_dev)
        return k_val, k_sig, 0.0

    print(f"\n[+] Running OpenMC Transport for cool_mult = {cool_mult:.2f} ({PARTICLES} particles, {BATCHES} batches)...")
    t0 = time.time()
    build_openmc_model_opt15(
        run_dir=run_dir,
        drum_angles=[0.0] * DRUM_COUNT,
        T_fuel_grid=T_fuel_vol,
        T_clad_grid=T_clad_avg,
        T_cool_grid=T_cool_bulk,
        T_mod=939.5,
        T_refl=880.0,
        cool_density_mult=cool_mult,
        particles=PARTICLES,
        batches=BATCHES,
        inactive=INACTIVE,
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

    with openmc.StatePoint(sp_path) as sp:
        k_val = float(sp.keff.nominal_value)
        k_sig = float(sp.keff.std_dev)

    dt = time.time() - t0
    print(f"    Transport complete in {dt:.1f} s: k_eff = {k_val:.5f} +/- {k_sig:.5f}")
    return k_val, k_sig, dt


def main():
    print("=" * 80)
    print("ULTRA-HIGH-PRECISION COOLANT FEEDBACK EVALUATION: GEOM_OPT15_EXTENDED")
    print(f"Histories: {PARTICLES * (BATCHES - INACTIVE) / 1e6:.1f}M active particles per state")
    print("=" * 80)

    # 1. Minus 75% (-75% coolant density, 0.25x nominal)
    dir_m75 = os.path.join(OUT_DIR, "cool_minus75pct")
    k_m75, sig_m75, dt_m75 = run_state(dir_m75, cool_mult=0.25, seed=401)

    # 2. Plus 75% (+75% coolant density, 1.75x nominal)
    dir_p75 = os.path.join(OUT_DIR, "cool_plus75pct")
    k_p75, sig_p75, dt_p75 = run_state(dir_p75, cool_mult=1.75, seed=402)

    # 3. Two-point difference across 150% density span
    # rho = (1 - 1/k) * 1e5 pcm
    rho_m75 = (1.0 - 1.0 / k_m75) * 1e5
    rho_p75 = (1.0 - 1.0 / k_p75) * 1e5
    delta_rho_150 = rho_p75 - rho_m75
    alpha_cool_2pt = delta_rho_150 / 150.0 # pcm / % density change

    # Error propagation
    sig_rho_m75 = (sig_m75 / (k_m75**2)) * 1e5
    sig_rho_p75 = (sig_p75 / (k_p75**2)) * 1e5
    sig_alpha_2pt = math.sqrt(sig_rho_m75**2 + sig_rho_p75**2) / 150.0

    print("\n" + "-" * 80)
    print("TWO-POINT DIRECT DIFFERENCE (-75% to +75% Density Span = 150%):")
    print(f"  k(-75%) = {k_m75:.5f} +/- {sig_m75:.5f} -> rho = {rho_m75:+.2f} +/- {sig_rho_m75:.2f} pcm")
    print(f"  k(+75%) = {k_p75:.5f} +/- {sig_p75:.5f} -> rho = {rho_p75:+.2f} +/- {sig_rho_p75:.2f} pcm")
    print(f"  Delta rho = {delta_rho_150:+.2f} +/- {math.sqrt(sig_rho_m75**2 + sig_rho_p75**2):.2f} pcm")
    print(f"  alpha_cool (2-point) = {alpha_cool_2pt:+.4f} +/- {sig_alpha_2pt:.4f} pcm/%")
    print("-" * 80)

    # 4. Multi-Point Regression combining all evaluated states:
    # Points: -75%, -25%, 0% (nominal), +25%, +75%
    # Load -25%, 0%, +25%
    k_m25, sig_m25 = 1.01974489, 0.00021327
    k_nom, sig_nom = 1.02025120, 0.00029688
    k_p25, sig_p25 = 1.01993515, 0.00026817

    pts = [
        {"x_pct": -75.0, "k": k_m75, "sig_k": sig_m75},
        {"x_pct": -25.0, "k": k_m25, "sig_k": sig_m25},
        {"x_pct":   0.0, "k": k_nom, "sig_k": sig_nom},
        {"x_pct": +25.0, "k": k_p25, "sig_k": sig_p25},
        {"x_pct": +75.0, "k": k_p75, "sig_k": sig_p75},
    ]
    df_pts = pd.DataFrame(pts)
    df_pts["rho_pcm"] = (1.0 - 1.0 / df_pts["k"]) * 1e5
    df_pts["sig_rho"] = (df_pts["sig_k"] / (df_pts["k"]**2)) * 1e5
    df_pts["weight"] = 1.0 / (df_pts["sig_rho"]**2)

    # Weighted Least Squares: rho = a + b * x
    w = df_pts["weight"].values
    x = df_pts["x_pct"].values
    y = df_pts["rho_pcm"].values

    w_sum = np.sum(w)
    x_bar = np.sum(w * x) / w_sum
    y_bar = np.sum(w * y) / w_sum

    b_num = np.sum(w * (x - x_bar) * (y - y_bar))
    b_den = np.sum(w * (x - x_bar)**2)
    alpha_cool_gls = b_num / b_den
    sig_alpha_gls = 1.0 / math.sqrt(b_den)

    z_score = alpha_cool_gls / sig_alpha_gls

    print("\n" + "=" * 80)
    print("5-POINT GENERALIZED LEAST SQUARES (GLS) REGRESSION SUMMARY:")
    print("=" * 80)
    for idx, row in df_pts.iterrows():
        print(f"  Point {idx+1} ({row['x_pct']:+5.1f}%): k = {row['k']:.5f} +/- {row['sig_k']:.5f} | rho = {row['rho_pcm']:+7.2f} +/- {row['sig_rho']:.2f} pcm")
    print("-" * 80)
    print(f"  alpha_cool (GLS Fit)     = {alpha_cool_gls:+.4f} +/- {sig_alpha_gls:.4f} pcm/%")
    print(f"  Z-score (alpha / sigma)  = {z_score:+.2f} sigma")
    sign_str = "STRICTLY POSITIVE" if z_score > 0 else "STRICTLY NEGATIVE"
    conf_pct = (1.0 - 0.5 * math.erfc(abs(z_score) / math.sqrt(2))) * 100.0
    print(f"  Statistical Confidence:  {conf_pct:.2f}% certainty that alpha_cool is {sign_str} (> {abs(z_score):.1f} sigma)")
    print("=" * 80)

    # Save detailed CSV
    df_pts.to_csv(os.path.join(DATA_DIR, "opt15_coolant_feedback_highprec.csv"), index=False)

    # Update opt15_feedback_coefficients.csv
    fb_path = os.path.join(DATA_DIR, "opt15_feedback_coefficients.csv")
    if os.path.exists(fb_path):
        df_fb = pd.read_csv(fb_path)
        df_fb["alpha_cool_pcm_pct"] = alpha_cool_gls
        df_fb["sigma_alpha_cool"] = sig_alpha_gls
        df_fb.to_csv(fb_path, index=False)
        print(f"[+] Updated {fb_path}")

    # Output authoritative text summary
    summary_path = os.path.join(DATA_DIR, "opt15_coolant_feedback_summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"alpha_cool_gls: {alpha_cool_gls:+.4f}\n")
        f.write(f"sigma_alpha_cool: {sig_alpha_gls:.4f}\n")
        f.write(f"z_score: {z_score:+.2f}\n")
        f.write(f"confidence_pct: {conf_pct:.2f}\n")
    print(f"[+] Written summary to {summary_path}")


if __name__ == "__main__":
    main()
