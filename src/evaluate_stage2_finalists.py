#!/usr/bin/env python3
"""
===============================================================================
STAGE 2 INTERMEDIATE EVALUATION & DOWN-SELECTION ENGINE: CANDIDATE OPT-10
===============================================================================
Takes the top finalists from Stage 1 screening and performs intermediate-fidelity
evaluation across all safety, neutronics, and mass dimensions:
  1. Operating Multiplication: k_out (BOL, 0 deg)
  2. All-in Shutdown Margin: SDM (180 deg) >= 2000 pcm
  3. Single Stuck-Drum Margin: SDM_stuck >= 1000 pcm
  4. Prompt Fuel Doppler: alpha_fuel (300 K to 1500 K) < 0 at 3-sigma
  5. Moderator Expansion: alpha_mod (300 K to 900 K) < 0 at 3-sigma
  6. Coolant Void Reactivity: Delta-rho_void < 0
  7. Exact Mass & Volume Accounting
  8. Multi-Criteria Pareto Scoring
===============================================================================
"""

import os
import sys
import csv
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import openmc

os.environ["PATH"] = "/home/aleri/miniforge3/envs/openmc-env/bin:" + os.environ.get("PATH", "")
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
os.environ["OMP_NUM_THREADS"] = "8"
openmc.config['cross_sections'] = "/home/aleri/wmp_data/cross_sections_wmp.xml"

# Import builders from screen_architectures
from screen_architectures import (
    ArchitectureCandidate,
    build_candidate_materials,
    build_candidate_geometry,
    compute_candidate_inventory
)


def run_stage2_evaluation(finalist_candidates, output_dir="results/stage2_evaluation"):
    os.makedirs(output_dir, exist_ok=True)
    work_dir = os.path.join(output_dir, "workspace")
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()

    # Moderate-high fidelity settings
    particles_std = 35000
    batches_std = 40
    inactive_std = 10

    # Higher particles for void and feedback
    particles_prec = 50000
    batches_prec = 50
    inactive_prec = 12

    results = []
    os.chdir(work_dir)

    for cand in finalist_candidates:
        print("\n" + "="*80)
        print(f"STAGE 2 DEEP EVALUATION: [{cand.cand_id}] {cand.name} ({cand.family})")
        print("="*80)
        t_start = time.time()

        sett_std = openmc.Settings()
        sett_std.run_mode = "eigenvalue"
        sett_std.particles = particles_std
        sett_std.batches = batches_std
        sett_std.inactive = inactive_std
        sett_std.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}
        sett_std.source = openmc.IndependentSource(space=openmc.stats.Point((0, 0, 0)))

        sett_prec = openmc.Settings()
        sett_prec.run_mode = "eigenvalue"
        sett_prec.particles = particles_prec
        sett_prec.batches = batches_prec
        sett_prec.inactive = inactive_prec
        sett_prec.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}
        sett_prec.source = openmc.IndependentSource(space=openmc.stats.Point((0, 0, 0)))

        # Inventory
        openmc.reset_auto_ids()
        mats_dummy, f_d, cl_d, co_d, mo_d, be_d, b_d = build_candidate_materials(cand)
        geom_dummy, pins, _ = build_candidate_geometry(cand, f_d, cl_d, co_d, mo_d, be_d, b_d)
        inv = compute_candidate_inventory(cand, pins)
        print(f"  Configuration: Pins={inv['n_pins']} | Mass={inv['total_mass_kg']:.2f} kg | Vol={inv['total_volume_L']:.2f} L | H/R={cand.aspect_ratio_H_R:.3f}")

        # 1. BOL Operating State (0 deg)
        print("  1. Operating State (Drums 0 deg) ...")
        openmc.reset_auto_ids()
        mats_out, f_o, cl_o, co_o, mo_o, be_o, b_o = build_candidate_materials(cand)
        geom_out, _, _ = build_candidate_geometry(cand, f_o, cl_o, co_o, mo_o, be_o, b_o, drum_angles=[0.0]*cand.number_of_drums)
        mod_out = openmc.Model(materials=mats_out, geometry=geom_out, settings=sett_std)
        mod_out.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{batches_std}.h5") as sp:
            k_out = float(sp.keff.nominal_value)
            s_out = float(sp.keff.std_dev)

        # 2. All-in Shutdown State (180 deg)
        print("  2. All-In Shutdown State (Drums 180 deg) ...")
        openmc.reset_auto_ids()
        mats_in, f_i, cl_i, co_i, mo_i, be_i, b_i = build_candidate_materials(cand)
        geom_in, _, _ = build_candidate_geometry(cand, f_i, cl_i, co_i, mo_i, be_i, b_i, drum_angles=[180.0]*cand.number_of_drums)
        mod_in = openmc.Model(materials=mats_in, geometry=geom_in, settings=sett_std)
        mod_in.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{batches_std}.h5") as sp:
            k_in = float(sp.keff.nominal_value)
            s_in = float(sp.keff.std_dev)

        sdm = (1.0 - k_in) * 1e5
        drum_worth = (1.0 / k_in - 1.0 / k_out) * 1e5

        # 3. 1-Stuck Drum State (Drum 0 at 0 deg, remaining at 180 deg)
        print("  3. 1-Stuck Drum Safety State ...")
        stuck_angles = [0.0] + [180.0] * (cand.number_of_drums - 1)
        openmc.reset_auto_ids()
        mats_stk, f_s, cl_s, co_s, mo_s, be_s, b_s = build_candidate_materials(cand)
        geom_stk, _, _ = build_candidate_geometry(cand, f_s, cl_s, co_s, mo_s, be_s, b_s, drum_angles=stuck_angles)
        mod_stk = openmc.Model(materials=mats_stk, geometry=geom_stk, settings=sett_std)
        mod_stk.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{batches_std}.h5") as sp:
            k_stuck = float(sp.keff.nominal_value)
            s_stuck = float(sp.keff.std_dev)

        sdm_stuck = (1.0 - k_stuck) * 1e5

        # 4. Prompt Fuel Doppler Feedback (300 K vs 1500 K)
        print("  4. Fuel Doppler Reactivity Feedback (300 K vs 1500 K) ...")
        # Cold Fuel 300 K
        openmc.reset_auto_ids()
        mats_fc, f_fc, cl_fc, co_fc, mo_fc, be_fc, b_fc = build_candidate_materials(cand)
        f_fc.temperature = 300.0
        geom_fc, _, _ = build_candidate_geometry(cand, f_fc, cl_fc, co_fc, mo_fc, be_fc, b_fc, drum_angles=[0.0]*cand.number_of_drums)
        mod_fc = openmc.Model(materials=mats_fc, geometry=geom_fc, settings=sett_std)
        mod_fc.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{batches_std}.h5") as sp:
            k_fc = float(sp.keff.nominal_value)
            s_fc = float(sp.keff.std_dev)

        # Hot Fuel 1500 K
        openmc.reset_auto_ids()
        mats_fh, f_fh, cl_fh, co_fh, mo_fh, be_fh, b_fh = build_candidate_materials(cand)
        f_fh.temperature = 1500.0
        geom_fh, _, _ = build_candidate_geometry(cand, f_fh, cl_fh, co_fh, mo_fh, be_fh, b_fh, drum_angles=[0.0]*cand.number_of_drums)
        mod_fh = openmc.Model(materials=mats_fh, geometry=geom_fh, settings=sett_std)
        mod_fh.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{batches_std}.h5") as sp:
            k_fh = float(sp.keff.nominal_value)
            s_fh = float(sp.keff.std_dev)

        alpha_fuel = (((k_fh - 1.0)/k_fh - (k_fc - 1.0)/k_fc) * 1e5) / 1200.0
        u_alpha_fuel = np.sqrt((s_fc/k_fc**2)**2 + (s_fh/k_fh**2)**2) / 1200.0 * 1e5

        # 5. Moderator Thermal Expansion Feedback (300 K vs 900 K)
        print("  5. Moderator Thermal Expansion Feedback (300 K vs 900 K) ...")
        # Cold 300 K
        rho_m3 = cand.moderator_density * 1.0108
        openmc.reset_auto_ids()
        mats_mc, f_mc, cl_mc, co_mc, mo_mc, be_mc, b_mc = build_candidate_materials(cand)
        mo_mc.set_density("g/cm3", rho_m3)
        mo_mc.temperature = 300.0
        geom_mc, _, _ = build_candidate_geometry(cand, f_mc, cl_mc, co_mc, mo_mc, be_mc, b_mc, drum_angles=[0.0]*cand.number_of_drums)
        mod_mc = openmc.Model(materials=mats_mc, geometry=geom_mc, settings=sett_std)
        mod_mc.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{batches_std}.h5") as sp:
            k_mc = float(sp.keff.nominal_value)
            s_mc = float(sp.keff.std_dev)

        # Hot 900 K
        rho_m9 = cand.moderator_density * 0.9892
        openmc.reset_auto_ids()
        mats_mh, f_mh, cl_mh, co_mh, mo_mh, be_mh, b_mh = build_candidate_materials(cand)
        mo_mh.set_density("g/cm3", rho_m9)
        mo_mh.temperature = 900.0
        geom_mh, _, _ = build_candidate_geometry(cand, f_mh, cl_mh, co_mh, mo_mh, be_mh, b_mh, drum_angles=[0.0]*cand.number_of_drums)
        mod_mh = openmc.Model(materials=mats_mh, geometry=geom_mh, settings=sett_std)
        mod_mh.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{batches_std}.h5") as sp:
            k_mh = float(sp.keff.nominal_value)
            s_mh = float(sp.keff.std_dev)

        alpha_mod = (((k_mh - 1.0)/k_mh - (k_mc - 1.0)/k_mc) * 1e5) / 600.0
        u_alpha_mod = np.sqrt((s_mc/k_mc**2)**2 + (s_mh/k_mh**2)**2) / 600.0 * 1e5

        # Total feedback (neglecting small coolant gas mass)
        alpha_tot = alpha_fuel + alpha_mod
        u_alpha_tot = np.sqrt(u_alpha_fuel**2 + u_alpha_mod**2)

        # 6. Coolant Void Reactivity (Precision Run)
        print("  6. Coolant Void Reactivity (Precision Run) ...")
        # Operating
        openmc.reset_auto_ids()
        mats_vo, f_vo, cl_vo, co_vo, mo_vo, be_vo, b_vo = build_candidate_materials(cand)
        geom_vo, _, _ = build_candidate_geometry(cand, f_vo, cl_vo, co_vo, mo_vo, be_vo, b_vo, drum_angles=[0.0]*cand.number_of_drums)
        mod_vo = openmc.Model(materials=mats_vo, geometry=geom_vo, settings=sett_prec)
        mod_vo.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{batches_prec}.h5") as sp:
            k_vo = float(sp.keff.nominal_value)
            s_vo = float(sp.keff.std_dev)

        # Voided
        openmc.reset_auto_ids()
        mats_vv, f_vv, cl_vv, co_vv, mo_vv, be_vv, b_vv = build_candidate_materials(cand, cool_density=1.0e-7)
        geom_vv, _, _ = build_candidate_geometry(cand, f_vv, cl_vv, co_vv, mo_vv, be_vv, b_vv, drum_angles=[0.0]*cand.number_of_drums)
        mod_vv = openmc.Model(materials=mats_vv, geometry=geom_vv, settings=sett_prec)
        mod_vv.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{batches_prec}.h5") as sp:
            k_vv = float(sp.keff.nominal_value)
            s_vv = float(sp.keff.std_dev)

        drho_void = (1.0 / k_vo - 1.0 / k_vv) * 1e5
        u_void = np.sqrt((s_vo/k_vo**2)**2 + (s_vv/k_vv**2)**2) * 1e5

        dt = time.time() - t_start
        print(f"  Summary [{cand.cand_id}]:")
        print(f"    k_out = {k_out:.5f} +/- {s_out:.5f} | k_in = {k_in:.5f} | k_stuck = {k_stuck:.5f}")
        print(f"    SDM = {sdm:.0f} pcm | SDM_stuck = {sdm_stuck:.0f} pcm | Worth = {drum_worth:.0f} pcm")
        print(f"    alpha_fuel = {alpha_fuel:+.2f} +/- {u_alpha_fuel:.2f} pcm/K (3-sig: {alpha_fuel+3*u_alpha_fuel:+.2f})")
        print(f"    alpha_mod  = {alpha_mod:+.2f} +/- {u_alpha_mod:.2f} pcm/K (3-sig: {alpha_mod+3*u_alpha_mod:+.2f})")
        print(f"    alpha_tot  = {alpha_tot:+.2f} +/- {u_alpha_tot:.2f} pcm/K (3-sig: {alpha_tot+3*u_alpha_tot:+.2f})")
        print(f"    drho_void  = {drho_void:+.1f} +/- {u_void:.1f} pcm (3-sig bound: {drho_void+3*u_void:+.1f})")
        print(f"    Mass = {inv['total_mass_kg']:.2f} kg | Volume = {inv['total_volume_L']:.2f} L | Evaluation Time = {dt:.1f}s")

        eval_res = {
            "cand_id": cand.cand_id,
            "family": cand.family,
            "name": cand.name,
            "mass_total_kg": inv["total_mass_kg"],
            "volume_total_L": inv["total_volume_L"],
            "pins": inv["n_pins"],
            "v_mod_v_fuel": inv["v_mod_v_fuel_ratio"],
            "k_out": k_out, "s_out": s_out,
            "k_in": k_in, "s_in": s_in,
            "k_stuck": k_stuck, "s_stuck": s_stuck,
            "sdm_pcm": sdm,
            "sdm_stuck_pcm": sdm_stuck,
            "drum_worth_pcm": drum_worth,
            "alpha_fuel_pcm_K": alpha_fuel, "u_alpha_fuel": u_alpha_fuel,
            "alpha_mod_pcm_K": alpha_mod, "u_alpha_mod": u_alpha_mod,
            "alpha_tot_pcm_K": alpha_tot, "u_alpha_tot": u_alpha_tot,
            "drho_void_pcm": drho_void, "u_void_pcm": u_void,
            "void_3sig_upper_pcm": drho_void + 3.0 * u_void,
            "alpha_tot_3sig_upper": alpha_tot + 3.0 * u_alpha_tot,
            "m_U235_kg": inv["m_U235_kg"],
            "m_fuel_kg": inv["m_fuel_kg"],
            "m_refl_kg": inv["m_refl_kg"],
            "eval_time_s": dt
        }
        results.append(eval_res)

    os.chdir(orig_dir)

    # Save to CSV
    csv_path = os.path.join(output_dir, "stage2_finalists_evaluation.csv")
    keys = list(results[0].keys())
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print(f"\nSaved Stage 2 results to {csv_path}")
    return results

if __name__ == "__main__":
    pass
