#!/usr/bin/env python3
"""
run_simulation_fourth_pipeline.py

Master pipeline script for Simulation 4 Campaign:
Runs both 4A (Base 4-Layer) and 4B (5-Layer with Al2O3 Liner) end-to-end,
then executes the multi-campaign comparison synthesis and figures generation.
"""

import os
import sys
import subprocess
import time

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/final_calculation/simulation_fourth"
PYTHON_EXEC = "/home/aleri/miniforge3/envs/openmc-env/bin/python"

steps = [
    # --- SIMULATION 4A (BASE 4-LAYER) ---
    ("4A: Coupled Picard Iteration", "sim4a_base/scripts/run_coupled_picard_4a.py", "sim4a_base"),
    ("4A: BOL Neutronics Suite", "sim4a_base/scripts/run_bol_neutronics_4a.py", "sim4a_base"),
    ("4A: Reactivity Feedbacks Suite", "sim4a_base/scripts/run_feedbacks_4a.py", "sim4a_base"),
    ("4A: Multi-Year Depletion Analysis", "sim4a_base/scripts/run_depletion_4a.py", "sim4a_base"),

    # --- SIMULATION 4B (5-LAYER WITH AL2O3 LINER) ---
    ("4B: Coupled Picard Iteration", "sim4b_al2o3/scripts/run_coupled_picard_4b.py", "sim4b_al2o3"),
    ("4B: BOL Neutronics Suite", "sim4b_al2o3/scripts/run_bol_neutronics_4b.py", "sim4b_al2o3"),
    ("4B: Reactivity Feedbacks Suite", "sim4b_al2o3/scripts/run_feedbacks_4b.py", "sim4b_al2o3"),
    ("4B: Multi-Year Depletion Analysis", "sim4b_al2o3/scripts/run_depletion_4b.py", "sim4b_al2o3"),

    # --- SYNTHESIS & COMPARISON ---
    ("Comparison: Data Synthesis", "comparison/scripts/run_full_comparison_study.py", "comparison"),
    ("Comparison: Figures Generation", "comparison/scripts/generate_comparison_figures_and_report.py", "comparison"),
    ("Comparison: Text Document Export", "comparison/scripts/generate_final_comparison_results_text.py", "comparison")
]

print("=" * 80)
print("STARTING SIMULATION 4 DUAL-CONFIGURATION PIPELINE (4A & 4B)")
print("=" * 80)

total_t0 = time.time()

for step_name, script_rel, cwd_rel in steps:
    script_path = os.path.join(BASE_DIR, script_rel)
    cwd_path = os.path.join(BASE_DIR, cwd_rel)
    print(f"\n{'='*80}")
    print(f"[RUNNING] {step_name} ({script_rel})")
    print(f"{'='*80}")
    t0 = time.time()

    cmd = [PYTHON_EXEC, script_path]
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"

    proc = subprocess.run(cmd, env=env, cwd=cwd_path)
    dt = time.time() - t0

    if proc.returncode != 0:
        print(f"\n[!] ERROR: Step '{step_name}' failed with returncode {proc.returncode}")
        sys.exit(proc.returncode)
    else:
        print(f"\n[+] SUCCESS: Step '{step_name}' completed in {dt:.1f} s")

total_dt = time.time() - total_t0
print("\n" + "=" * 80)
print(f"ALL SIMULATION 4 STEPS COMPLETED SUCCESSFULLY IN {total_dt/60:.1f} MINUTES!")
print("=" * 80)
