"""
High-Statistics Coolant Void Reactivity Evaluation for Candidate Opt-6
======================================================================
Evaluates the coolant void reactivity at 25,000 particles/batch and 60 batches
to resolve the statistical error below 40 pcm and verify confidence intervals.
"""

import os
import sys
import numpy as np
import openmc

os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
if "/home/aleri/miniforge3/envs/openmc-env/bin" not in os.environ.get("PATH", ""):
    os.environ["PATH"] = f"/home/aleri/miniforge3/envs/openmc-env/bin:{os.environ.get('PATH', '')}"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import ReactorParameters, build_materials, build_geometry, build_settings

def run_void_evaluation(work_dir="scratch/test_void_opt6"):
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    p_op = ReactorParameters(
        core_radius=17.90, active_core_height=35.00,
        fuel_radius=0.58, cladding_thickness=0.05,
        coolant_channel_thickness=0.045, moderator_fuel_pitch=1.615,
        radial_reflector_thickness=9.60, axial_reflector_thickness=7.50,
        control_drum_radius=4.30, absorber_thickness=1.45,
        absorber_sector_angular_width=140.0, b10_enrichment=0.90,
        coolant_state='operating',
        particles=25000, batches=60, inactive=12
    )

    print("\n--- 1. Evaluating Operating Coolant State (25k particles, 60 batches) ---")
    openmc.reset_auto_ids()
    mats_op, f_op, cl_op, co_op, mo_op, be_op, b_op = build_materials(p_op)
    g_op, _, _ = build_geometry(p_op, f_op, cl_op, co_op, mo_op, be_op, b_op)
    mod_op = openmc.Model(materials=mats_op, geometry=g_op, settings=build_settings(p_op))
    mod_op.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p_op.batches}.h5") as sp:
        k_op = float(sp.keff.nominal_value)
        s_op = float(sp.keff.std_dev)
    print(f"  Operating keff = {k_op:.5f} +/- {s_op:.5f}")

    print("\n--- 2. Evaluating Voided Coolant State (25k particles, 60 batches) ---")
    p_vd = ReactorParameters(
        core_radius=17.90, active_core_height=35.00,
        fuel_radius=0.58, cladding_thickness=0.05,
        coolant_channel_thickness=0.045, moderator_fuel_pitch=1.615,
        radial_reflector_thickness=9.60, axial_reflector_thickness=7.50,
        control_drum_radius=4.30, absorber_thickness=1.45,
        absorber_sector_angular_width=140.0, b10_enrichment=0.90,
        coolant_state='voided',
        particles=25000, batches=60, inactive=12
    )
    p_vd.coolant_density = 1.0e-7
    openmc.reset_auto_ids()
    mats_vd, f_vd, cl_vd, co_vd, mo_vd, be_vd, b_vd = build_materials(p_vd)
    g_vd, _, _ = build_geometry(p_vd, f_vd, cl_vd, co_vd, mo_vd, be_vd, b_vd)
    mod_vd = openmc.Model(materials=mats_vd, geometry=g_vd, settings=build_settings(p_vd))
    mod_vd.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p_vd.batches}.h5") as sp:
        k_vd = float(sp.keff.nominal_value)
        s_vd = float(sp.keff.std_dev)
    print(f"  Voided keff    = {k_vd:.5f} +/- {s_vd:.5f}")

    drho_void = ((k_vd - 1.0)/k_vd - (k_op - 1.0)/k_op) * 1e5
    u_void = np.sqrt((s_vd/k_vd**2)**2 + (s_op/k_op**2)**2) * 1e5

    print("\n" + "=" * 80)
    print("COOLANT VOID REACTIVITY SUMMARY:")
    print("=" * 80)
    print(f"  Delta-rho_void = {drho_void:+.1f} +/- {u_void:.1f} pcm (1-sigma)")
    print(f"  2-sigma interval: [{drho_void - 2*u_void:+.1f}, {drho_void + 2*u_void:+.1f}] pcm")
    print(f"  3-sigma interval: [{drho_void - 3*u_void:+.1f}, {drho_void + 3*u_void:+.1f}] pcm")

    os.chdir(orig_dir)
    return drho_void, u_void

if __name__ == "__main__":
    run_void_evaluation()
