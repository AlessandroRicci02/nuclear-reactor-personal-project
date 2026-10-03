"""
Candidate Screening Runner for Opt-5 Ultra-Compact 400 kWth Reactor
===================================================================
Evaluates neutronic performance and exact mass/volume across candidate
geometries to select the optimal Opt-5 design point:
- Core radius sweep: 15.0 cm to 16.5 cm
- Active core height: 42.0 cm to 44.0 cm
- Control drum radius & absorber sizing
- 400 kWth power rating requirements:
    * k_out(BOL) >= 1.035 - 1.045
    * k_in(BOL) <= 0.970 (SDM >= 3000 pcm)
    * k_stuck <= 0.985 (SDM_stk >= 1500 pcm)
    * Total mass minimized (< 550 kg target)
"""

import os
import sys
import csv
import numpy as np
import openmc

os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
if "/home/aleri/miniforge3/envs/openmc-env/bin" not in os.environ.get("PATH", ""):
    os.environ["PATH"] = f"/home/aleri/miniforge3/envs/openmc-env/bin:{os.environ.get('PATH', '')}"

from model import ReactorParameters, build_materials, build_geometry, build_settings, calculate_inventory

def count_pins(p):
    grid_dim = int(np.ceil(2 * p.core_radius / p.moderator_fuel_pitch))
    if grid_dim % 2 == 0:
        grid_dim += 1
    count = 0
    for i in range(grid_dim):
        for j in range(grid_dim):
            x = -grid_dim * p.moderator_fuel_pitch / 2 + (j + 0.5) * p.moderator_fuel_pitch
            y = -grid_dim * p.moderator_fuel_pitch / 2 + (i + 0.5) * p.moderator_fuel_pitch
            if x**2 + y**2 <= (p.core_radius - p.channel_outer_radius)**2:
                count += 1
    return count

def evaluate_candidate(cand, work_dir="results/opt5/workspace"):
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)
    try:
        p = ReactorParameters(
            core_radius=cand["r_core"],
            active_core_height=cand["h_core"],
            moderator_fuel_pitch=cand.get("pitch", 1.635),
            radial_reflector_thickness=cand.get("t_refl", 9.60),
            axial_reflector_thickness=cand.get("t_axial", 7.60),
            control_drum_radius=cand.get("r_drum", 4.30),
            number_of_control_drums=8,
            drum_clearance=0.20,
            absorber_thickness=1.45,
            absorber_sector_angular_width=140.0,
            b10_enrichment=0.90,
            particles=5000,
            batches=28,
            inactive=8
        )
        pins = count_pins(p)
        inv = calculate_inventory(p, pins)

        # 1. k_out (drums out, 0°)
        openmc.reset_auto_ids()
        p.drum_angle_deg = 0.0
        mats, f, cl, co, m, be, b = build_materials(p)
        geom, _, _ = build_geometry(p, f, cl, co, m, be, b)
        sett = build_settings(p)
        mod_out = openmc.Model(materials=mats, geometry=geom, settings=sett)
        mod_out.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
            k_out = float(sp.keff.nominal_value)
            s_out = float(sp.keff.std_dev)

        # 2. k_in (drums in, 180°)
        openmc.reset_auto_ids()
        p.drum_angle_deg = 180.0
        mats_in, f_in, cl_in, co_in, m_in, be_in, b_in = build_materials(p)
        geom_in, _, _ = build_geometry(p, f_in, cl_in, co_in, m_in, be_in, b_in)
        sett_in = build_settings(p)
        mod_in = openmc.Model(materials=mats_in, geometry=geom_in, settings=sett_in)
        mod_in.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
            k_in = float(sp.keff.nominal_value)
            s_in = float(sp.keff.std_dev)

        # 3. k_stuck (1 stuck at 0°, 7 at 180°)
        openmc.reset_auto_ids()
        stuck_angles = [0.0] + [180.0] * 7
        mats_stk, f_stk, cl_stk, co_stk, m_stk, be_stk, b_stk = build_materials(p)
        geom_stk, _, _ = build_geometry(p, f_stk, cl_stk, co_stk, m_stk, be_stk, b_stk, drum_angles=stuck_angles)
        sett_stk = build_settings(p)
        mod_stk = openmc.Model(materials=mats_stk, geometry=geom_stk, settings=sett_stk)
        mod_stk.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
            k_stk = float(sp.keff.nominal_value)
            s_stk = float(sp.keff.std_dev)

        drum_worth_pcm = (1.0 / k_in - 1.0 / k_out) * 1e5
        sdm_normal_pcm = (1.0 - k_in) / k_in * 1e5
        sdm_stuck_pcm = (1.0 - k_stk) / k_stk * 1e5

        res = {
            "name": cand["name"],
            "r_core": cand["r_core"],
            "h_core": cand["h_core"],
            "pins": pins,
            "fuel_mass_kg": inv["total_fuel_mass_kg"],
            "total_mass_kg": inv["total_mass_kg"],
            "total_volume_L": inv["total_volume_L"],
            "delta_mass_kg": inv["total_mass_kg"] - 622.68,
            "k_out": k_out, "s_out": s_out,
            "k_in": k_in, "s_in": s_in,
            "k_stuck": k_stk, "s_stuck": s_stk,
            "drum_worth_pcm": drum_worth_pcm,
            "sdm_normal_pcm": sdm_normal_pcm,
            "sdm_stuck_pcm": sdm_stuck_pcm
        }
        return res
    finally:
        os.chdir(orig_dir)

def main():
    print("=" * 80)
    print("CANDIDATE SCREENING FOR OPT-5 ULTRA-COMPACT 400 kWth SPACE REACTOR")
    print("=" * 80)

    candidates = [
        {"name": "Cand-5A (R=16.5, H=44.0)", "r_core": 16.50, "h_core": 44.00, "pitch": 1.635, "t_refl": 9.80, "t_axial": 7.80, "r_drum": 4.40},
        {"name": "Cand-5B (R=16.0, H=43.5)", "r_core": 16.00, "h_core": 43.50, "pitch": 1.635, "t_refl": 9.60, "t_axial": 7.60, "r_drum": 4.30},
        {"name": "Cand-5C (R=15.5, H=43.0)", "r_core": 15.50, "h_core": 43.00, "pitch": 1.635, "t_refl": 9.50, "t_axial": 7.50, "r_drum": 4.20},
        {"name": "Cand-5D (R=15.0, H=42.0)", "r_core": 15.00, "h_core": 42.00, "pitch": 1.635, "t_refl": 9.50, "t_axial": 7.50, "r_drum": 4.10},
    ]

    results = []
    for cand in candidates:
        print(f"\nEvaluating {cand['name']}...")
        r = evaluate_candidate(cand)
        results.append(r)
        print(f"  Pins: {r['pins']} | Total Mass: {r['total_mass_kg']:.2f} kg ({r['delta_mass_kg']:+.2f} kg vs Opt-4)")
        print(f"  k_out: {r['k_out']:.5f} +/- {r['s_out']:.5f} (Excess = {(r['k_out']-1)/r['k_out']*1e5:+.0f} pcm)")
        print(f"  k_in : {r['k_in']:.5f} +/- {r['s_in']:.5f} (SDM = {r['sdm_normal_pcm']:4.0f} pcm)")
        print(f"  k_stk: {r['k_stuck']:.5f} +/- {r['s_stuck']:.5f} (SDM_stk = {r['sdm_stuck_pcm']:4.0f} pcm)")
        print(f"  Drum Worth: {r['drum_worth_pcm']:.0f} pcm")

    # Export CSV summary
    csv_path = "results/opt5/candidate_screening_summary.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Candidate", "R_core [cm]", "H_core [cm]", "Pins", "Fuel Mass [kg]",
            "Total Mass [kg]", "Delta Mass [kg]", "Volume [L]",
            "k_out", "s_out", "k_in", "s_in", "k_stuck", "s_stuck",
            "Drum Worth [pcm]", "Normal SDM [pcm]", "Stuck SDM [pcm]"
        ])
        for r in results:
            writer.writerow([
                r["name"], f"{r['r_core']:.2f}", f"{r['h_core']:.2f}", r["pins"],
                f"{r['fuel_mass_kg']:.2f}", f"{r['total_mass_kg']:.2f}", f"{r['delta_mass_kg']:+.2f}",
                f"{r['total_volume_L']:.2f}", f"{r['k_out']:.5f}", f"{r['s_out']:.5f}",
                f"{r['k_in']:.5f}", f"{r['s_in']:.5f}", f"{r['k_stuck']:.5f}", f"{r['s_stuck']:.5f}",
                f"{r['drum_worth_pcm']:.0f}", f"{r['sdm_normal_pcm']:.0f}", f"{r['sdm_stuck_pcm']:.0f}"
            ])
    print(f"\nSaved screening summary to: {csv_path}")

if __name__ == "__main__":
    main()
