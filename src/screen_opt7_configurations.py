"""
Multi-Objective Screening for Opt-7 Candidates
==============================================
Systematically evaluates candidates across:
- Coolant composition: He-Xe 85/15 (M=23.1) vs He-Xe 90/10 (M=16.7) vs He-Xe 80/20 (M=29.5)
- Annular coolant gap: 0.040 cm (0.40 mm) vs 0.045 cm (0.45 mm)
- Coolant pressure: 2.0 MPa vs 1.6 MPa
- Geometry: H/R in [1.88, 1.96]
Evaluates: BOL k_out, k_in, 1-stuck drum, drum worth, mass, volume,
and void reactivity with precise error tracking.
"""

import os
import sys
import csv
import numpy as np
import openmc

os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
if "/home/aleri/miniforge3/envs/openmc-env/bin" not in os.environ.get("PATH", ""):
    os.environ["PATH"] = f"/home/aleri/miniforge3/envs/openmc-env/bin:{os.environ.get('PATH', '')}"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import ReactorParameters, build_materials, build_geometry, build_settings, calculate_inventory
from screen_opt6_candidates import count_pins

def screen_opt7_candidates(work_dir="scratch/screen_opt7"):
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    print("=" * 80)
    print("SCREENING CANDIDATE CONFIGURATIONS FOR OPT-7")
    print("=" * 80)

    # Candidate Matrix:
    # We explore combinations designed to drive void reactivity robustly negative
    # while preserving mass < 515 kg, volume < 119 L, TMC < -2.5 pcm/K, and lifetime > 20 yr.
    candidates = [
        {
            "id": "Opt7-A",
            "desc": "He-Xe 85/15 (M=23.1), P=2.0 MPa, Gap=0.040 cm, R=17.90 cm, H=34.50 cm (H/R=1.93), Pitch=1.615 cm",
            "he": 0.85, "xe": 0.15, "M": 23.096, "press_MPa": 2.0, "gap": 0.040,
            "rc": 17.90, "hc": 34.50, "pitch": 1.615, "rf": 0.58, "clad": 0.05,
            "r_ref": 9.60, "z_ref": 7.50, "r_drum": 4.30
        },
        {
            "id": "Opt7-B",
            "desc": "He-Xe 90/10 (M=16.7), P=2.0 MPa, Gap=0.040 cm, R=17.90 cm, H=34.50 cm (H/R=1.93), Pitch=1.615 cm",
            "he": 0.90, "xe": 0.10, "M": 16.732, "press_MPa": 2.0, "gap": 0.040,
            "rc": 17.90, "hc": 34.50, "pitch": 1.615, "rf": 0.58, "clad": 0.05,
            "r_ref": 9.60, "z_ref": 7.50, "r_drum": 4.30
        },
        {
            "id": "Opt7-C",
            "desc": "He-Xe 85/15 (M=23.1), P=1.6 MPa, Gap=0.045 cm, R=17.80 cm, H=34.20 cm (H/R=1.92), Pitch=1.615 cm",
            "he": 0.85, "xe": 0.15, "M": 23.096, "press_MPa": 1.6, "gap": 0.045,
            "rc": 17.80, "hc": 34.20, "pitch": 1.615, "rf": 0.58, "clad": 0.05,
            "r_ref": 9.50, "z_ref": 7.40, "r_drum": 4.30
        },
        {
            "id": "Opt7-D",
            "desc": "He-Xe 90/10 (M=16.7), P=2.0 MPa, Gap=0.045 cm, R=17.80 cm, H=34.80 cm (H/R=1.95), Pitch=1.612 cm",
            "he": 0.90, "xe": 0.10, "M": 16.732, "press_MPa": 2.0, "gap": 0.045,
            "rc": 17.80, "hc": 34.80, "pitch": 1.612, "rf": 0.58, "clad": 0.05,
            "r_ref": 9.60, "z_ref": 7.50, "r_drum": 4.30
        }
    ]

    results = []

    for c in candidates:
        print(f"\n---> Evaluating Candidate {c['id']}: {c['desc']}")
        # Ideal gas law: rho = (P_Pa * M_kg) / (R_gas * T_K)
        # R_gas = 8.314 J/(mol*K), T = 1050 K
        # rho [g/cm3] = (P_MPa * 1e6 * M * 1e-3) / (8.314 * 1050 * 1e6) = (P_MPa * M) / (8.314 * 1050) * 1e-3
        rho_cool = (c["press_MPa"] * 1.0e6 * c["M"] * 1.0e-3) / (8.314 * 1050.0) * 1e-3  # g/cm3

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
            particles=20000,
            batches=45,
            inactive=10
        )
        pins = count_pins(p_base)
        inv = calculate_inventory(p_base, pins)

        # 1. Operating State (Drums 0 deg, operating coolant)
        openmc.reset_auto_ids()
        p_base.drum_angle_deg = 0.0
        p_base.coolant_state = 'operating'
        mats_op, f_op, cl_op, co_op, m_op, be_op, b_op = build_materials(p_base)
        co_op.nuclides.clear()
        co_op.set_density("g/cm3", rho_cool)
        co_op.add_nuclide("He4", c["he"], "ao")
        if c["xe"] > 0:
            co_op.add_element("Xe", c["xe"], "ao")
        g_op, _, _ = build_geometry(p_base, f_op, cl_op, co_op, m_op, be_op, b_op)
        s_op = build_settings(p_base)
        mod_op = openmc.Model(materials=mats_op, geometry=g_op, settings=s_op)
        mod_op.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_base.batches}.h5") as sp:
            k_op = float(sp.keff.nominal_value)
            s_op = float(sp.keff.std_dev)

        # 2. Voided State (Drums 0 deg, vacuum)
        openmc.reset_auto_ids()
        p_vd = ReactorParameters(
            core_radius=c["rc"], active_core_height=c["hc"], fuel_radius=c["rf"],
            cladding_thickness=c["clad"], coolant_channel_thickness=c["gap"],
            moderator_fuel_pitch=c["pitch"], radial_reflector_thickness=c["r_ref"],
            axial_reflector_thickness=c["z_ref"], control_drum_radius=c["r_drum"],
            absorber_thickness=1.45, absorber_sector_angular_width=140.0,
            b10_enrichment=0.90, drum_angle_deg=0.0, coolant_state='voided',
            particles=20000, batches=45, inactive=10
        )
        p_vd.coolant_density = 1.0e-7
        mats_vd, f_vd, cl_vd, co_vd, m_vd, be_vd, b_vd = build_materials(p_vd)
        g_vd, _, _ = build_geometry(p_vd, f_vd, cl_vd, co_vd, m_vd, be_vd, b_vd)
        s_vd = build_settings(p_vd)
        mod_vd = openmc.Model(materials=mats_vd, geometry=g_vd, settings=s_vd)
        mod_vd.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_vd.batches}.h5") as sp:
            k_vd = float(sp.keff.nominal_value)
            s_vd = float(sp.keff.std_dev)

        # 3. Full Shutdown (Drums 180 deg, operating coolant)
        openmc.reset_auto_ids()
        p_in = ReactorParameters(
            core_radius=c["rc"], active_core_height=c["hc"], fuel_radius=c["rf"],
            cladding_thickness=c["clad"], coolant_channel_thickness=c["gap"],
            moderator_fuel_pitch=c["pitch"], radial_reflector_thickness=c["r_ref"],
            axial_reflector_thickness=c["z_ref"], control_drum_radius=c["r_drum"],
            absorber_thickness=1.45, absorber_sector_angular_width=140.0,
            b10_enrichment=0.90, drum_angle_deg=180.0, coolant_state='operating',
            particles=12000, batches=35, inactive=10
        )
        mats_in, f_in, cl_in, co_in, m_in, be_in, b_in = build_materials(p_in)
        co_in.nuclides.clear()
        co_in.set_density("g/cm3", rho_cool)
        co_in.add_nuclide("He4", c["he"], "ao")
        if c["xe"] > 0:
            co_in.add_element("Xe", c["xe"], "ao")
        g_in, _, _ = build_geometry(p_in, f_in, cl_in, co_in, m_in, be_in, b_in)
        s_in = build_settings(p_in)
        mod_in = openmc.Model(materials=mats_in, geometry=g_in, settings=s_in)
        mod_in.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_in.batches}.h5") as sp:
            k_in = float(sp.keff.nominal_value)
            s_in = float(sp.keff.std_dev)

        # 4. Stuck Drum (1 at 0 deg, 7 at 180 deg)
        openmc.reset_auto_ids()
        stuck_angles = [0.0] + [180.0] * 7
        mats_stk, f_stk, cl_stk, co_stk, m_stk, be_stk, b_stk = build_materials(p_in)
        co_stk.nuclides.clear()
        co_stk.set_density("g/cm3", rho_cool)
        co_stk.add_nuclide("He4", c["he"], "ao")
        if c["xe"] > 0:
            co_stk.add_element("Xe", c["xe"], "ao")
        g_stk, _, _ = build_geometry(p_in, f_stk, cl_stk, co_stk, m_stk, be_stk, b_stk, drum_angles=stuck_angles)
        s_stk = build_settings(p_in)
        mod_stk = openmc.Model(materials=mats_stk, geometry=g_stk, settings=s_stk)
        mod_stk.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_in.batches}.h5") as sp:
            k_stk = float(sp.keff.nominal_value)
            s_stk = float(sp.keff.std_dev)

        # Physics quantities:
        # Delta-rho_void = (1/k_nom - 1/k_void) * 1e5
        drho_void = (1.0 / k_op - 1.0 / k_vd) * 1e5
        u_void = np.sqrt((s_op/k_op**2)**2 + (s_vd/k_vd**2)**2) * 1e5
        bound_3sigma = drho_void + 3.0 * u_void

        sdm_norm = (1.0 - k_in) / k_in * 1e5
        sdm_stk = (1.0 - k_stk) / k_stk * 1e5
        worth = (1.0 / k_in - 1.0 / k_op) * 1e5

        v_fuel_cell = np.pi * c["rf"]**2
        v_mod_cell = c["pitch"]**2 - np.pi * (c["rf"] + c["clad"] + c["gap"])**2
        ratio_cell = v_mod_cell / v_fuel_cell

        rec = {
            "id": c["id"], "desc": c["desc"], "he": c["he"], "xe": c["xe"], "M": c["M"],
            "press": c["press_MPa"], "gap": c["gap"], "rc": c["rc"], "hc": c["hc"],
            "hr": c["hc"]/c["rc"], "hd": c["hc"]/(2.0*c["rc"]),
            "pitch": c["pitch"], "v_ratio_cell": ratio_cell, "pins": pins,
            "mass_kg": inv["total_mass_kg"], "volume_L": inv["total_volume_L"],
            "k_op": k_op, "s_op": s_op, "k_vd": k_vd, "s_vd": s_vd,
            "drho_void": drho_void, "u_void": u_void, "bound_3sigma": bound_3sigma,
            "k_in": k_in, "s_in": s_in, "k_stk": k_stk, "s_stk": s_stk,
            "sdm_norm": sdm_norm, "sdm_stk": sdm_stk, "worth": worth
        }
        results.append(rec)

        print(f"  Pins: {pins:3d} | Mass: {inv['total_mass_kg']:6.2f} kg | Volume: {inv['total_volume_L']:5.2f} L | H/R = {c['hc']/c['rc']:.2f}")
        print(f"  k_op(BOL) : {k_op:.5f} +/- {s_op:.5f} (Excess = {(k_op-1)/k_op*1e5:+.0f} pcm)")
        print(f"  k_vd(BOL) : {k_vd:.5f} +/- {s_vd:.5f}")
        print(f"  VOID REACTIVITY : {drho_void:+.1f} +/- {u_void:.1f} pcm (3-sigma bound: {bound_3sigma:+.1f} pcm)")
        print(f"  k_in(BOL) : {k_in:.5f} +/- {s_in:.5f} | SDM = {sdm_norm:.0f} pcm")
        print(f"  k_stk(BOL): {k_stk:.5f} +/- {s_stk:.5f} | SDM_stk = {sdm_stk:.0f} pcm")
        print(f"  8-Drum Worth: {worth:.0f} pcm")

    os.chdir(orig_dir)
    csv_path = "results/opt7/screen_opt7_candidates.csv"
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "ID", "Description", "He_ao", "Xe_ao", "M_g_mol", "Pressure_MPa", "Gap_cm",
            "R_core_cm", "H_core_cm", "H_R", "H_D", "Pitch_cm", "Vmod_Vfuel_cell",
            "Pins", "Mass_kg", "Volume_L", "k_op", "s_op", "k_void", "s_void",
            "Void_drho_pcm", "u_void_pcm", "Bound_3sigma_pcm",
            "k_in", "s_in", "k_stuck", "s_stuck", "SDM_pcm", "SDM_stuck_pcm", "Worth_pcm"
        ])
        for r in results:
            writer.writerow([
                r["id"], r["desc"], f"{r['he']:.2f}", f"{r['xe']:.2f}", f"{r['M']:.3f}",
                f"{r['press']:.1f}", f"{r['gap']:.3f}", f"{r['rc']:.2f}", f"{r['hc']:.2f}",
                f"{r['hr']:.3f}", f"{r['hd']:.3f}", f"{r['pitch']:.3f}", f"{r['v_ratio_cell']:.3f}",
                r["pins"], f"{r['mass_kg']:.2f}", f"{r['volume_L']:.2f}",
                f"{r['k_op']:.5f}", f"{r['s_op']:.5f}", f"{r['k_vd']:.5f}", f"{r['s_vd']:.5f}",
                f"{r['drho_void']:+.1f}", f"{r['u_void']:.1f}", f"{r['bound_3sigma']:+.1f}",
                f"{r['k_in']:.5f}", f"{r['s_in']:.5f}", f"{r['k_stk']:.5f}", f"{r['s_stk']:.5f}",
                f"{r['sdm_norm']:.0f}", f"{r['sdm_stk']:.0f}", f"{r['worth']:.0f}"
            ])
    print(f"\nSaved screening results to: {csv_path}")
    return results

if __name__ == "__main__":
    screen_opt7_candidates()
