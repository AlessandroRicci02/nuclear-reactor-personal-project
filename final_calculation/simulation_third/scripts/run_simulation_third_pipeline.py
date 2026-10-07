#!/usr/bin/env python3
"""
run_simulation_third_pipeline.py

Master pipeline script for Third Simulation Campaign:
Executes sequentially:
  1. BOL Neutronics (high-statistics)
  2. Reactivity Feedback Coefficients (high-statistics)
  3. Multi-Year Depletion Analysis (2.5-year steps to 25 years)
  4. Thermal-Hydraulic & Brayton Systems Analysis
  5. 300 DPI Publication-Quality Figures Generation
  6. Authoritative Final Results Document Export
"""

import os
import sys
import subprocess
import time

BASE_DIR = "/home/aleri/OpenMC_projects/first_test/final_calculation/simulation_third"
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")
PYTHON_EXEC = "/home/aleri/miniforge3/envs/openmc-env/bin/python"

scripts_to_run = [
    ("BOL Neutronics Suite", "run_bol_neutronics_900_1150k.py"),
    ("Reactivity Feedbacks Suite", "run_feedbacks_900_1150k.py"),
    ("Multi-Year Depletion Analysis", "run_depletion_2.5yr_900_1150k.py"),
    ("Thermal-Hydraulics & Systems Analysis", "export_th_and_thermodynamics_third.py"),
    ("Figures Generation (300 DPI)", "generate_figures_and_report_third.py"),
    ("Final Results Document Export", "generate_final_results_text.py")
]

print("=" * 80)
print("STARTING SIMULATION 3 AUTOMATED EXECUTION PIPELINE")
print("=" * 80)

for step_name, script_name in scripts_to_run:
    script_path = os.path.join(SCRIPTS_DIR, script_name)
    print(f"\n{'='*80}")
    print(f"[RUNNING STEP] {step_name} ({script_name})")
    print(f"{'='*80}")
    t0 = time.time()
    
    cmd = [PYTHON_EXEC, script_path]
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    
    proc = subprocess.run(cmd, env=env, cwd=BASE_DIR)
    dt = time.time() - t0
    
    if proc.returncode != 0:
        print(f"\n[!] ERROR: Step '{step_name}' failed with returncode {proc.returncode}")
        sys.exit(proc.returncode)
    else:
        print(f"\n[+] SUCCESS: Step '{step_name}' completed in {dt:.1f} s")

print("\n" + "=" * 80)
print("SIMULATION 3 PIPELINE COMPLETED SUCCESSFULLY!")
print("=" * 80)
