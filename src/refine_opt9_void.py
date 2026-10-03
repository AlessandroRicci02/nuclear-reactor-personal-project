#!/usr/bin/env python3
"""
Refined Ultra-High Statistics Void Reactivity Run for Opt-9:
Particles: 500,000
Batches: 130 (20 inactive, 110 active)
Target: sigma <= 16 pcm so that Delta-rho_void + 3*sigma < 0
"""

import os
import sys
import csv
import numpy as np
import openmc

os.environ["PATH"] = "/home/aleri/miniforge3/envs/openmc-env/bin:" + os.environ.get("PATH", "")
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
os.environ["OMP_NUM_THREADS"] = "8"
openmc.config['cross_sections'] = "/home/aleri/wmp_data/cross_sections_wmp.xml"

from run_opt9_master import Opt9Parameters, build_opt9_materials, build_opt9_geometry, build_opt9_settings

def run_refined_void():
    p = Opt9Parameters()
    work_dir = "results/opt9/workspace_void_refine"
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()

    try:
        os.chdir(work_dir)
        settings = build_opt9_settings(p, particles=500000, batches=130, inactive=20)

        # 1. Operating State (Pure He 2.5 MPa)
        print("Running Operating State (500k particles x 130 batches) ...")
        openmc.reset_auto_ids()
        m_op, f_op, cl_op, co_op, mo_op, be_op, b_op = build_opt9_materials(p)
        g_op, _, _ = build_opt9_geometry(p, f_op, cl_op, co_op, mo_op, be_op, b_op, drum_angles=[0.0]*8)
        mod_op = openmc.Model(materials=m_op, geometry=g_op, settings=settings)
        mod_op.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{settings.batches}.h5") as sp:
            k_op = float(sp.keff.nominal_value)
            s_op = float(sp.keff.std_dev)

        # 2. Voided State (Vacuum)
        print("Running Voided State (500k particles x 130 batches) ...")
        openmc.reset_auto_ids()
        m_vd, f_vd, cl_vd, co_vd, mo_vd, be_vd, b_vd = build_opt9_materials(p, cool_density=1.0e-7)
        g_vd, _, _ = build_opt9_geometry(p, f_vd, cl_vd, co_vd, mo_vd, be_vd, b_vd, drum_angles=[0.0]*8)
        mod_vd = openmc.Model(materials=m_vd, geometry=g_vd, settings=settings)
        mod_vd.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{settings.batches}.h5") as sp:
            k_vd = float(sp.keff.nominal_value)
            s_vd = float(sp.keff.std_dev)

        drho_void = (1.0 / k_op - 1.0 / k_vd) * 1e5
        u_void = np.sqrt((s_op/k_op**2)**2 + (s_vd/k_vd**2)**2) * 1e5
        bound_3sigma = drho_void + 3.0 * u_void

        print(f"\n=======================================================")
        print(f"REFINED OPT-9 VOID RESULTS:")
        print(f"  Operating keff = {k_op:.5f} +/- {s_op:.5f}")
        print(f"  Voided keff    = {k_vd:.5f} +/- {s_vd:.5f}")
        print(f"  Delta-rho_void = {drho_void:+.1f} +/- {u_void:.1f} pcm")
        print(f"  3-sigma bound  = {bound_3sigma:+.1f} pcm")
        print(f"=======================================================\n")

        # Update void CSV
        void_csv = os.path.join(orig_dir, "results/opt9/void/void_reactivity.csv")
        with open(void_csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["State", "Coolant", "Density_g_cm3", "keff", "std_dev", "Delta_rho_pcm", "u_pcm", "Bound_3sigma_pcm"])
            writer.writerow(["Operating", "Pure He (2.5 MPa)", f"{p.coolant_density:.6f}", f"{k_op:.5f}", f"{s_op:.5f}", "0.0", "0.0", "0.0"])
            writer.writerow(["Voided (Vacuum)", "Void / Axial Ducts", "0.000000", f"{k_vd:.5f}", f"{s_vd:.5f}", f"{drho_void:+.1f}", f"{u_void:.1f}", f"{bound_3sigma:+.1f}"])

        return k_op, s_op, k_vd, s_vd, drho_void, u_void, bound_3sigma

    finally:
        os.chdir(orig_dir)

if __name__ == "__main__":
    run_refined_void()
