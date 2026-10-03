#!/usr/bin/env python3
"""
Multi-objective ranking and downselection of Stage 1 screening candidates.
"""
import os
import csv
import numpy as np

def score_candidates(csv_path="results/architecture_screening.csv"):
    if not os.path.exists(csv_path):
        print(f"{csv_path} not found.")
        return None

    candidates = []
    with open(csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            candidates.append(row)

    print("\n" + "="*110)
    print("STAGE 1 CANDIDATE ARCHITECTURE SCREENING RESULTS & COMPREHENSIVE COMPARISON")
    print("="*110)
    print(f"{'ID':<3} | {'Family':<24} | {'Name':<26} | {'Mass':<6} | {'Vol':<5} | {'k_out':<8} | {'k_in':<8} | {'SDM':<6} | {'Worth':<6} | {'Void (pcm)':<14}")
    print("-" * 110)

    for c in candidates:
        m = float(c["mass_total_kg"])
        v = float(c["volume_total_L"])
        ko = float(c["k_out"])
        ki = float(c["k_in"])
        sdm = float(c["sdm_pcm"])
        w = float(c["drum_worth_pcm"])
        vd = float(c["drho_void_pcm"])
        uvd = float(c["u_void_pcm"])
        print(f"{c['cand_id']:<3} | {c['family'][:24]:<24} | {c['name'][:26]:<26} | {m:<6.1f} | {v:<5.1f} | {ko:<8.5f} | {ki:<8.5f} | {sdm:<6.0f} | {w:<6.0f} | {vd:+5.1f} +/- {uvd:<4.1f}")

    print("="*110)

if __name__ == "__main__":
    score_candidates()
