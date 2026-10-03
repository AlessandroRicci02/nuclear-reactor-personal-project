#!/usr/bin/env python3
"""
Test Candidate Opt-9 configurations using the exact validated OpenMC CSG geometry engine.
"""

import os
import sys
import numpy as np
import openmc

os.environ["PATH"] = "/home/aleri/miniforge3/envs/openmc-env/bin:" + os.environ.get("PATH", "")
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
os.environ["OMP_NUM_THREADS"] = "8"
openmc.config['cross_sections'] = "/home/aleri/wmp_data/cross_sections_wmp.xml"

from run_opt8_master import Opt8Parameters, build_opt8_materials, build_opt8_geometry, build_opt8_settings

def evaluate_candidate(r_core, h_core, r_refl_rad=10.2, h_refl_ax=9.5):
    p = Opt8Parameters()
    p.core_radius = r_core
    p.active_core_height = h_core
    p.radial_reflector_thickness = r_refl_rad
    p.axial_reflector_thickness = h_refl_ax
    p.radial_reflector_outer_radius = p.core_radius + p.radial_reflector_thickness
    p.total_reactor_height = p.active_core_height + 2 * p.axial_reflector_thickness
    p.drum_radial_position = p.core_radius + p.drum_clearance + p.control_drum_radius

    # Mass and volume
    vol_tot = np.pi * p.radial_reflector_outer_radius**2 * p.total_reactor_height * 1e-3
    v_core = np.pi * p.core_radius**2 * p.active_core_height
    v_refl_rad = np.pi * (p.radial_reflector_outer_radius**2 - p.core_radius**2) * p.total_reactor_height
    v_refl_ax = np.pi * p.core_radius**2 * 2 * p.axial_reflector_thickness
    # Opt-8 calibrated mass density formula
    mass_est = 628.49 * (vol_tot / 154.06)

    sett = build_opt8_settings(p, particles=40000, batches=45, inactive=10)

    # 1. Operating Drums OUT (0 deg)
    mats1, f1, cl1, co1, mo1, re1, ab1 = build_opt8_materials(p)
    geom1, pins1, _ = build_opt8_geometry(p, f1, cl1, co1, mo1, re1, ab1, drum_angles=[0.0]*8)
    model1 = openmc.Model(materials=mats1, geometry=geom1, settings=sett)
    model1.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
        k_out = float(sp.keff.nominal_value)
        s_out = float(sp.keff.std_dev)

    # 2. Operating Drums OUT - VOIDED
    mats2, f2, cl2, co2, mo2, re2, ab2 = build_opt8_materials(p, cool_density=1.0e-7)
    geom2, pins2, _ = build_opt8_geometry(p, f2, cl2, co2, mo2, re2, ab2, drum_angles=[0.0]*8)
    model2 = openmc.Model(materials=mats2, geometry=geom2, settings=sett)
    model2.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
        k_vd = float(sp.keff.nominal_value)
        s_vd = float(sp.keff.std_dev)

    # 3. Operating Drums IN (180 deg)
    mats3, f3, cl3, co3, mo3, re3, ab3 = build_opt8_materials(p)
    geom3, pins3, _ = build_opt8_geometry(p, f3, cl3, co3, mo3, re3, ab3, drum_angles=[180.0]*8)
    model3 = openmc.Model(materials=mats3, geometry=geom3, settings=sett)
    model3.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
        k_in = float(sp.keff.nominal_value)
        s_in = float(sp.keff.std_dev)

    drho_void = (1.0/k_out - 1.0/k_vd) * 1e5
    u_void = np.sqrt((s_out/k_out**2)**2 + (s_vd/k_vd**2)**2) * 1e5
    bound_3s_void = drho_void + 3.0 * u_void
    sdm = (1.0 - k_in) * 1e5
    worth = (1.0/k_in - 1.0/k_out) * 1e5

    print(f"\n======================================================================")
    print(f"CANDIDATE: R={r_core:.2f} cm, H={h_core:.2f} cm (H/R={h_core/r_core:.2f})")
    print(f"  Pins={len(pins1)} | Vol={vol_tot:.1f} L | Mass_est={mass_est:.1f} kg")
    print(f"  k_out (0°)   = {k_out:.5f} +/- {s_out:.5f}")
    print(f"  k_in  (180°) = {k_in:.5f} +/- {s_in:.5f}")
    print(f"  SDM          = {sdm:.0f} pcm | Drum Worth = {worth:.0f} pcm")
    print(f"  Void drho    = {drho_void:+.1f} +/- {u_void:.1f} pcm | 3-sigma bound = {bound_3s_void:+.1f} pcm")
    print(f"======================================================================\n")

if __name__ == "__main__":
    # Test 3 carefully selected candidates:
    # 1. R=18.0 cm (365 pins), H=38.0 cm (H/R=2.11)
    evaluate_candidate(18.00, 38.00)
    # 2. R=18.2 cm (373 pins), H=37.5 cm (H/R=2.06)
    evaluate_candidate(18.20, 37.50)
    # 3. R=18.2 cm (373 pins), H=36.0 cm (H/R=1.98)
    evaluate_candidate(18.20, 36.00)
