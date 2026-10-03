"""
Deep Statistical Physics & Temperature Feedback Analysis for Opt-7 Candidates
=============================================================================
Evaluates Candidate Opt-7 (R=17.80 cm, H=34.80 cm, H/R=1.955, P=1.612 cm, 349 pins):
1. Very high statistics void reactivity across:
   - He-Xe 90/10 (M = 16.732 g/mol, standard high-conductivity Brayton gas)
   - He-Xe 95/5  (M = 10.367 g/mol)
   - Pure Helium (M = 4.003 g/mol)
   Using 60,000 particles/batch and 70 batches (3.6M active histories) to push sigma_drho <= 35 pcm!
2. Full temperature reactivity feedback (Doppler, MTC, Coolant, and total TMC)
   to ensure unconditional stability (alpha_total + 3*sigma < 0).
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

def run_deep_void_analysis(work_dir="scratch/deep_void_opt7"):
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    print("=" * 80)
    print("HIGH-STATISTICS VOID REACTIVITY (60,000 PARTICLES/BATCH, 70 BATCHES)")
    print("=" * 80)

    p_geom = ReactorParameters(
        core_radius=17.80, active_core_height=34.80,
        fuel_radius=0.58, cladding_thickness=0.05,
        coolant_channel_thickness=0.045, moderator_fuel_pitch=1.612,
        radial_reflector_thickness=9.60, axial_reflector_thickness=7.50,
        control_drum_radius=4.30, absorber_thickness=1.45,
        absorber_sector_angular_width=140.0, b10_enrichment=0.90,
        particles=60000, batches=70, inactive=12
    )

    # 1. Vacuum Voided State
    print("\n--- 1. Simulating Voided State (3.6M active histories) ---")
    openmc.reset_auto_ids()
    p_vd = ReactorParameters(
        core_radius=p_geom.core_radius, active_core_height=p_geom.active_core_height,
        fuel_radius=p_geom.fuel_radius, cladding_thickness=p_geom.cladding_thickness,
        coolant_channel_thickness=p_geom.coolant_channel_thickness,
        moderator_fuel_pitch=p_geom.moderator_fuel_pitch,
        radial_reflector_thickness=p_geom.radial_reflector_thickness,
        axial_reflector_thickness=p_geom.axial_reflector_thickness,
        control_drum_radius=p_geom.control_drum_radius,
        absorber_thickness=p_geom.absorber_thickness,
        absorber_sector_angular_width=p_geom.absorber_sector_angular_width,
        b10_enrichment=p_geom.b10_enrichment,
        coolant_state='voided',
        particles=60000, batches=70, inactive=12
    )
    p_vd.coolant_density = 1.0e-7
    mats_v, f_v, cl_v, co_v, mo_v, be_v, b_v = build_materials(p_vd)
    g_v, _, _ = build_geometry(p_vd, f_v, cl_v, co_v, mo_v, be_v, b_v)
    s_v = build_settings(p_vd)
    mod_v = openmc.Model(materials=mats_v, geometry=g_v, settings=s_v)
    mod_v.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p_vd.batches}.h5") as sp:
        k_void = float(sp.keff.nominal_value)
        s_void = float(sp.keff.std_dev)
    print(f"Voided Reference: k_void = {k_void:.5f} +/- {s_void:.5f}")

    gas_options = [
        {"name": "He-Xe 90/10 (M=16.7)", "he": 0.90, "xe": 0.10, "M": 16.732, "rho": 0.00383},
        {"name": "He-Xe 95/5  (M=10.4)", "he": 0.95, "xe": 0.05, "M": 10.367, "rho": 0.00238},
        {"name": "Pure He    (M=4.0)",  "he": 1.00, "xe": 0.00, "M":  4.003, "rho": 0.00092},
    ]

    void_results = []
    for g in gas_options:
        print(f"\n--- Simulating Operating State: {g['name']} ---")
        openmc.reset_auto_ids()
        mats, f, cl, co, m, be, b = build_materials(p_geom)
        co.nuclides.clear()
        co.set_density("g/cm3", g["rho"])
        co.add_nuclide("He4", g["he"], "ao")
        if g["xe"] > 0:
            co.add_element("Xe", g["xe"], "ao")
        g_op, _, _ = build_geometry(p_geom, f, cl, co, m, be, b)
        s_op = build_settings(p_geom)
        mod_op = openmc.Model(materials=mats, geometry=g_op, settings=s_op)
        mod_op.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_geom.batches}.h5") as sp:
            k_op = float(sp.keff.nominal_value)
            s_op = float(sp.keff.std_dev)

        # Delta-rho = (1/k_nom - 1/k_void) * 1e5
        drho = (1.0 / k_op - 1.0 / k_void) * 1e5
        u_drho = np.sqrt((s_op/k_op**2)**2 + (s_void/k_void**2)**2) * 1e5
        bound_3sigma = drho + 3.0 * u_drho

        rec = {
            "name": g["name"], "he": g["he"], "xe": g["xe"], "M": g["M"], "rho": g["rho"],
            "k_op": k_op, "s_op": s_op, "k_void": k_void, "s_void": s_void,
            "drho": drho, "u_drho": u_drho, "bound_3sigma": bound_3sigma
        }
        void_results.append(rec)
        print(f"  k_op = {k_op:.5f} +/- {s_op:.5f} | k_void = {k_void:.5f} +/- {s_void:.5f}")
        print(f"  Delta-rho_void = {drho:+.1f} +/- {u_drho:.1f} pcm")
        print(f"  1-sigma: [{drho - u_drho:+.1f}, {drho + u_drho:+.1f}] pcm")
        print(f"  2-sigma: [{drho - 2*u_drho:+.1f}, {drho + 2*u_drho:+.1f}] pcm")
        print(f"  3-sigma: [{drho - 3*u_drho:+.1f}, {drho + 3*u_drho:+.1f}] pcm")
        print(f"  3-sigma upper bound = {bound_3sigma:+.1f} pcm")

    os.chdir(orig_dir)
    return void_results

def run_opt7_temperature_feedback(gas_choice="He-Xe 90/10", he=0.90, xe=0.10, rho_cool=0.00383, work_dir="scratch/tmc_opt7"):
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    print("\n" + "=" * 80)
    print(f"OPT-7 CONTINUOUS-ENERGY TEMPERATURE REACTIVITY FEEDBACK ({gas_choice})")
    print("=" * 80)

    p = ReactorParameters(
        core_radius=17.80, active_core_height=34.80,
        fuel_radius=0.58, cladding_thickness=0.05,
        coolant_channel_thickness=0.045, moderator_fuel_pitch=1.612,
        radial_reflector_thickness=9.60, axial_reflector_thickness=7.50,
        control_drum_radius=4.30, absorber_thickness=1.45,
        absorber_sector_angular_width=140.0, b10_enrichment=0.90,
        particles=20000, batches=50, inactive=10
    )

    # 1. Fuel Doppler: 300 K vs 1500 K
    print("\n--- 1. Evaluating Fuel Doppler (300 K vs 1500 K) ---")
    openmc.reset_auto_ids()
    mats_f3, f3, cl, co3, mo, be, b = build_materials(p)
    f3.temperature = 300.0
    co3.nuclides.clear()
    co3.set_density("g/cm3", rho_cool)
    co3.add_nuclide("He4", he, "ao")
    if xe > 0:
        co3.add_element("Xe", xe, "ao")
    g_f3, _, _ = build_geometry(p, f3, cl, co3, mo, be, b)
    mod_f3 = openmc.Model(materials=mats_f3, geometry=g_f3, settings=build_settings(p))
    mod_f3.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_f3 = float(sp.keff.nominal_value)
        s_f3 = float(sp.keff.std_dev)

    openmc.reset_auto_ids()
    mats_f15, f15, cl, co15, mo, be, b = build_materials(p)
    f15.temperature = 1500.0
    co15.nuclides.clear()
    co15.set_density("g/cm3", rho_cool)
    co15.add_nuclide("He4", he, "ao")
    if xe > 0:
        co15.add_element("Xe", xe, "ao")
    g_f15, _, _ = build_geometry(p, f15, cl, co15, mo, be, b)
    mod_f15 = openmc.Model(materials=mats_f15, geometry=g_f15, settings=build_settings(p))
    mod_f15.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_f15 = float(sp.keff.nominal_value)
        s_f15 = float(sp.keff.std_dev)

    alpha_fuel = (((k_f15 - 1.0)/k_f15 - (k_f3 - 1.0)/k_f3) * 1e5) / 1200.0
    u_alpha_fuel = np.sqrt((s_f3/k_f3**2)**2 + (s_f15/k_f15**2)**2) / 1200.0 * 1e5
    print(f"  alpha_fuel = {alpha_fuel:+.2f} +/- {u_alpha_fuel:.2f} pcm/K (1-sigma)")

    # 2. Moderator Feedback: 300 K vs 900 K
    print("\n--- 2. Evaluating Moderator Feedback (300 K vs 900 K) ---")
    openmc.reset_auto_ids()
    mats_m3, f, cl, co_m3, m3, be, b = build_materials(p)
    m3.temperature = 300.0
    co_m3.nuclides.clear()
    co_m3.set_density("g/cm3", rho_cool)
    co_m3.add_nuclide("He4", he, "ao")
    if xe > 0:
        co_m3.add_element("Xe", xe, "ao")
    g_m3, _, _ = build_geometry(p, f, cl, co_m3, m3, be, b)
    mod_m3 = openmc.Model(materials=mats_m3, geometry=g_m3, settings=build_settings(p))
    mod_m3.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_m3 = float(sp.keff.nominal_value)
        s_m3 = float(sp.keff.std_dev)

    openmc.reset_auto_ids()
    mats_m9, f, cl, co_m9, m9, be, b = build_materials(p)
    m9.temperature = 900.0
    co_m9.nuclides.clear()
    co_m9.set_density("g/cm3", rho_cool)
    co_m9.add_nuclide("He4", he, "ao")
    if xe > 0:
        co_m9.add_element("Xe", xe, "ao")
    g_m9, _, _ = build_geometry(p, f, cl, co_m9, m9, be, b)
    mod_m9 = openmc.Model(materials=mats_m9, geometry=g_m9, settings=build_settings(p))
    mod_m9.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_m9 = float(sp.keff.nominal_value)
        s_m9 = float(sp.keff.std_dev)

    alpha_mod = (((k_m9 - 1.0)/k_m9 - (k_m3 - 1.0)/k_m3) * 1e5) / 600.0
    u_alpha_mod = np.sqrt((s_m3/k_m3**2)**2 + (s_m9/k_m9**2)**2) / 600.0 * 1e5
    print(f"  alpha_mod  = {alpha_mod:+.2f} +/- {u_alpha_mod:.2f} pcm/K (1-sigma)")

    # 3. Coolant Feedback: cold (900 K) vs hot (1300 K)
    print("\n--- 3. Evaluating Coolant Feedback (900 K vs 1300 K) ---")
    # At 2.0 MPa:
    # T = 900 K -> rho = (2e6 * M * 1e-3) / (8.314 * 900) * 1e-3 = rho_1050 * (1050/900) = 1.1667 * rho_1050
    # T = 1300 K -> rho = rho_1050 * (1050/1300) = 0.8077 * rho_1050
    rho_cold = rho_cool * (1050.0 / 900.0)
    rho_hot  = rho_cool * (1050.0 / 1300.0)

    openmc.reset_auto_ids()
    mats_cc, f_c, cl_c, co_c, mo_c, be_c, b_c = build_materials(p)
    co_c.nuclides.clear()
    co_c.set_density("g/cm3", rho_cold)
    co_c.add_nuclide("He4", he, "ao")
    if xe > 0:
        co_c.add_element("Xe", xe, "ao")
    g_cc, _, _ = build_geometry(p, f_c, cl_c, co_c, mo_c, be_c, b_c)
    mod_cc = openmc.Model(materials=mats_cc, geometry=g_cc, settings=build_settings(p))
    mod_cc.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_cc = float(sp.keff.nominal_value)
        s_cc = float(sp.keff.std_dev)

    openmc.reset_auto_ids()
    mats_ch, f_h, cl_h, co_h, mo_h, be_h, b_h = build_materials(p)
    co_h.nuclides.clear()
    co_h.set_density("g/cm3", rho_hot)
    co_h.add_nuclide("He4", he, "ao")
    if xe > 0:
        co_h.add_element("Xe", xe, "ao")
    g_ch, _, _ = build_geometry(p, f_h, cl_h, co_h, mo_h, be_h, b_h)
    mod_ch = openmc.Model(materials=mats_ch, geometry=g_ch, settings=build_settings(p))
    mod_ch.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_ch = float(sp.keff.nominal_value)
        s_ch = float(sp.keff.std_dev)

    alpha_cool = (((k_ch - 1.0)/k_ch - (k_cc - 1.0)/k_cc) * 1e5) / 400.0
    u_alpha_cool = np.sqrt((s_cc/k_cc**2)**2 + (s_ch/k_ch**2)**2) / 400.0 * 1e5
    print(f"  alpha_cool = {alpha_cool:+.2f} +/- {u_alpha_cool:.2f} pcm/K (1-sigma)")

    alpha_tot = alpha_fuel + alpha_mod + alpha_cool
    u_tot = np.sqrt(u_alpha_fuel**2 + u_alpha_mod**2 + u_alpha_cool**2)

    print("\n" + "=" * 80)
    print("TOTAL TEMPERATURE COEFFICIENT SUMMARY:")
    print("=" * 80)
    print(f"  alpha_fuel : {alpha_fuel:+.2f} +/- {u_alpha_fuel:.2f} pcm/K (1-sigma)")
    print(f"  alpha_mod  : {alpha_mod:+.2f} +/- {u_alpha_mod:.2f} pcm/K (1-sigma)")
    print(f"  alpha_cool : {alpha_cool:+.2f} +/- {u_alpha_cool:.2f} pcm/K (1-sigma)")
    print(f"  alpha_total: {alpha_tot:+.2f} +/- {u_tot:.2f} pcm/K (1-sigma)")
    print(f"  3-sigma upper bound for total TMC: {alpha_tot + 3*u_tot:+.2f} pcm/K")

    os.chdir(orig_dir)
    return alpha_fuel, u_alpha_fuel, alpha_mod, u_alpha_mod, alpha_cool, u_alpha_cool, alpha_tot, u_tot

if __name__ == "__main__":
    void_res = run_deep_void_analysis()
    tmc_res = run_opt7_temperature_feedback()
