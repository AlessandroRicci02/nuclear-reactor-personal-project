"""
Screening of Opt-6 NEP Reactor Candidates
=========================================
Systematically screens candidate configurations combining:
1. Optimal finite cylinder buckling ratio H/R in [1.80, 2.05]
2. Under-moderation V_mod/V_fuel in [0.85, 1.05] (Ade et al. 2022)
3. Void reactivity sign & statistical confidence (target: Delta-rho + 3*sigma < 0)
4. Cold shutdown margin (SDM > 2000 pcm) and stuck drum margin (SDM_stk > 1000 pcm)
5. Total assembly mass < 514 kg and volume < 120 L
"""

import os
import sys
import csv
import numpy as np
import openmc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
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

def screen_candidates(work_dir="scratch/screen_opt6"):
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    print("=" * 80)
    print("SCREENING CANDIDATE CONFIGURATIONS FOR OPT-6")
    print("=" * 80)

    # We evaluate 4 carefully chosen geometries around H/R = 1.85 - 2.0
    # and pitch = 1.56 - 1.60 cm (under-moderated)
    candidates = [
        {
            "id": "Opt6-Cand1",
            "desc": "H/R = 1.90, P = 1.58 cm, R = 17.50 cm, H = 33.25 cm (V_mod/V_fuel = 0.98)",
            "rc": 17.50, "hc": 33.25, "pitch": 1.58, "rf": 0.58, "clad": 0.05, "gap": 0.045,
            "r_ref": 9.40, "z_ref": 7.20, "r_drum": 4.20
        },
        {
            "id": "Opt6-Cand2",
            "desc": "H/R = 1.95, P = 1.57 cm, R = 17.40 cm, H = 33.90 cm (V_mod/V_fuel = 0.95)",
            "rc": 17.40, "hc": 33.90, "pitch": 1.57, "rf": 0.58, "clad": 0.05, "gap": 0.045,
            "r_ref": 9.30, "z_ref": 7.20, "r_drum": 4.20
        },
        {
            "id": "Opt6-Cand3",
            "desc": "H/R = 1.85, P = 1.56 cm, R = 17.80 cm, H = 32.90 cm (V_mod/V_fuel = 0.93)",
            "rc": 17.80, "hc": 32.90, "pitch": 1.56, "rf": 0.58, "clad": 0.05, "gap": 0.045,
            "r_ref": 9.30, "z_ref": 7.00, "r_drum": 4.20
        },
        {
            "id": "Opt6-Cand4",
            "desc": "H/R = 2.00, P = 1.59 cm, R = 17.20 cm, H = 34.40 cm (V_mod/V_fuel = 1.00)",
            "rc": 17.20, "hc": 34.40, "pitch": 1.59, "rf": 0.58, "clad": 0.05, "gap": 0.045,
            "r_ref": 9.40, "z_ref": 7.20, "r_drum": 4.20
        }
    ]

    results = []

    for c in candidates:
        print(f"\n---> Evaluating {c['id']}: {c['desc']}")

        # Parameter object
        p_base = ReactorParameters(
            core_radius=c["rc"],
            active_core_height=c["hc"],
            fuel_radius=c["rf"],
            cladding_thickness=c["clad"],
            coolant_channel_thickness=c["gap"],
            moderator_fuel_pitch=c["pitch"],
            radial_reflector_thickness=c["r_ref"],
            axial_reflector_thickness=c["z_ref"],
            control_drum_radius=c["r_drum"],
            absorber_thickness=1.45,
            absorber_sector_angular_width=140.0,
            b10_enrichment=0.90,
            particles=10000,
            batches=35,
            inactive=10
        )
        pins = count_pins(p_base)
        inv = calculate_inventory(p_base, pins)

        # 1. BOL Uncontrolled Out (0 deg)
        openmc.reset_auto_ids()
        p_base.drum_angle_deg = 0.0
        p_base.coolant_state = 'operating'
        mats_out, f, cl, co, m, be, b = build_materials(p_base)
        geom_out, _, _ = build_geometry(p_base, f, cl, co, m, be, b)
        sett_out = build_settings(p_base)
        mod_out = openmc.Model(materials=mats_out, geometry=geom_out, settings=sett_out)
        mod_out.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_base.batches}.h5") as sp:
            k_out = float(sp.keff.nominal_value)
            s_out = float(sp.keff.std_dev)

        # 2. BOL Full Shutdown (180 deg)
        openmc.reset_auto_ids()
        p_in = ReactorParameters(
            core_radius=c["rc"], active_core_height=c["hc"], fuel_radius=c["rf"],
            cladding_thickness=c["clad"], coolant_channel_thickness=c["gap"],
            moderator_fuel_pitch=c["pitch"], radial_reflector_thickness=c["r_ref"],
            axial_reflector_thickness=c["z_ref"], control_drum_radius=c["r_drum"],
            absorber_thickness=1.45, absorber_sector_angular_width=140.0,
            b10_enrichment=0.90, drum_angle_deg=180.0, coolant_state='operating',
            particles=10000, batches=35, inactive=10
        )
        mats_in, f_in, cl_in, co_in, m_in, be_in, b_in = build_materials(p_in)
        geom_in, _, _ = build_geometry(p_in, f_in, cl_in, co_in, m_in, be_in, b_in)
        sett_in = build_settings(p_in)
        mod_in = openmc.Model(materials=mats_in, geometry=geom_in, settings=sett_in)
        mod_in.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_in.batches}.h5") as sp:
            k_in = float(sp.keff.nominal_value)
            s_in = float(sp.keff.std_dev)

        # 3. BOL 1 Stuck Drum Out (1 drum at 0 deg, 7 at 180 deg)
        openmc.reset_auto_ids()
        stuck_angles = [0.0] + [180.0] * 7
        mats_stk, f_stk, cl_stk, co_stk, m_stk, be_stk, b_stk = build_materials(p_in)
        geom_stk, _, _ = build_geometry(p_in, f_stk, cl_stk, co_stk, m_stk, be_stk, b_stk, drum_angles=stuck_angles)
        sett_stk = build_settings(p_in)
        mod_stk = openmc.Model(materials=mats_stk, geometry=geom_stk, settings=sett_stk)
        mod_stk.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_in.batches}.h5") as sp:
            k_stuck = float(sp.keff.nominal_value)
            s_stuck = float(sp.keff.std_dev)

        # 4. Voided State (coolant_state = 'voided')
        openmc.reset_auto_ids()
        p_void = ReactorParameters(
            core_radius=c["rc"], active_core_height=c["hc"], fuel_radius=c["rf"],
            cladding_thickness=c["clad"], coolant_channel_thickness=c["gap"],
            moderator_fuel_pitch=c["pitch"], radial_reflector_thickness=c["r_ref"],
            axial_reflector_thickness=c["z_ref"], control_drum_radius=c["r_drum"],
            absorber_thickness=1.45, absorber_sector_angular_width=140.0,
            b10_enrichment=0.90, drum_angle_deg=0.0, coolant_state='voided',
            particles=10000, batches=35, inactive=10
        )
        mats_vd, f_vd, cl_vd, co_vd, m_vd, be_vd, b_vd = build_materials(p_void)
        geom_vd, _, _ = build_geometry(p_void, f_vd, cl_vd, co_vd, m_vd, be_vd, b_vd)
        sett_vd = build_settings(p_void)
        mod_vd = openmc.Model(materials=mats_vd, geometry=geom_vd, settings=sett_vd)
        mod_vd.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_void.batches}.h5") as sp:
            k_void = float(sp.keff.nominal_value)
            s_void = float(sp.keff.std_dev)

        # Physics derivations
        sdm_norm = (1.0 - k_in) / k_in * 1e5
        sdm_stk = (1.0 - k_stuck) / k_stuck * 1e5
        drum_worth = (1.0 / k_in - 1.0 / k_out) * 1e5

        drho_void = ((k_void - 1.0)/k_void - (k_out - 1.0)/k_out) * 1e5
        u_void = np.sqrt((s_void/k_void**2)**2 + (s_out/k_out**2)**2) * 1e5
        conf_3sigma_void = drho_void + 3.0 * u_void

        v_fuel_pin = np.pi * c["rf"]**2
        v_cell = c["pitch"]**2
        v_channel = np.pi * (c["rf"] + c["clad"] + c["gap"])**2
        v_mod_cell = v_cell - v_channel
        v_ratio = v_mod_cell / v_fuel_pin

        rec = {
            "id": c["id"],
            "desc": c["desc"],
            "v_ratio": v_ratio,
            "pins": pins,
            "mass_kg": inv["total_mass_kg"],
            "volume_L": inv["total_volume_L"],
            "k_out": k_out, "s_out": s_out,
            "k_in": k_in, "s_in": s_in,
            "k_stuck": k_stuck, "s_stuck": s_stuck,
            "sdm_norm_pcm": sdm_norm,
            "sdm_stk_pcm": sdm_stk,
            "drum_worth_pcm": drum_worth,
            "k_void": k_void, "s_void": s_void,
            "drho_void_pcm": drho_void,
            "u_void_pcm": u_void,
            "conf_3sigma_void_pcm": conf_3sigma_void
        }
        results.append(rec)

        print(f"  Pins: {pins:3d} | Mass: {inv['total_mass_kg']:6.2f} kg | Volume: {inv['total_volume_L']:5.2f} L")
        print(f"  k_out(BOL): {k_out:.5f} +/- {s_out:.5f} | Excess: {(k_out-1)/k_out*1e5:+.0f} pcm")
        print(f"  k_in(BOL) : {k_in:.5f} +/- {s_in:.5f} | SDM: {sdm_norm:5.0f} pcm")
        print(f"  k_stk(BOL): {k_stuck:.5f} +/- {s_stuck:.5f} | SDM_stk: {sdm_stk:5.0f} pcm")
        print(f"  8-Drum Worth: {drum_worth:5.0f} pcm ({drum_worth/1000:4.2f}% Delta-k/k)")
        print(f"  Void drho : {drho_void:+.1f} +/- {u_void:.1f} pcm | 3-sigma: {conf_3sigma_void:+.1f} pcm")

    os.chdir(orig_dir)
    csv_path = "results/opt6/screen_opt6_candidates.csv"
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Candidate", "Vmod/Vfuel", "Pins", "Mass [kg]", "Volume [L]",
            "k_out", "s_out", "k_in", "s_in", "k_stuck", "s_stuck",
            "SDM [pcm]", "SDM_stk [pcm]", "Drum Worth [pcm]",
            "k_void", "s_void", "Void drho [pcm]", "u_void [pcm]", "3-sigma upper bound [pcm]"
        ])
        for r in results:
            writer.writerow([
                r["id"], f"{r['v_ratio']:.3f}", r["pins"], f"{r['mass_kg']:.2f}", f"{r['volume_L']:.2f}",
                f"{r['k_out']:.5f}", f"{r['s_out']:.5f}", f"{r['k_in']:.5f}", f"{r['s_in']:.5f}",
                f"{r['k_stuck']:.5f}", f"{r['s_stuck']:.5f}",
                f"{r['sdm_norm_pcm']:.0f}", f"{r['sdm_stk_pcm']:.0f}", f"{r['drum_worth_pcm']:.0f}",
                f"{r['k_void']:.5f}", f"{r['s_void']:.5f}",
                f"{r['drho_void_pcm']:+.1f}", f"{r['u_void_pcm']:.1f}", f"{r['conf_3sigma_void_pcm']:+.1f}"
            ])
    print(f"\nSaved screening results to: {csv_path}")
    return results

if __name__ == "__main__":
    screen_candidates()
