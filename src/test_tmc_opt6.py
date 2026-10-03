"""
Test Temperature Coefficients of Reactivity for Candidate Opt-6
================================================================
Continuous-energy WMP evaluation of:
1. Fuel Doppler coefficient (300 K to 1500 K)
2. Moderator temperature coefficient (300 K to 900 K)
3. Coolant temperature coefficient (900 K to 1300 K)
4. Total temperature coefficient of reactivity with 1-sigma, 2-sigma, 3-sigma bounds
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

def run_tmc_test(work_dir="scratch/test_tmc_opt6"):
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    p = ReactorParameters(
        core_radius=17.90, active_core_height=35.00,
        fuel_radius=0.58, cladding_thickness=0.05,
        coolant_channel_thickness=0.045, moderator_fuel_pitch=1.615,
        radial_reflector_thickness=9.60, axial_reflector_thickness=7.50,
        control_drum_radius=4.30, absorber_thickness=1.45,
        absorber_sector_angular_width=140.0, b10_enrichment=0.90,
        particles=15000, batches=45, inactive=10
    )

    # 1. Fuel Doppler: 300 K vs 1500 K
    print("\n--- 1. Evaluating Fuel Doppler (300 K vs 1500 K) ---")
    openmc.reset_auto_ids()
    mats_f3, f3, cl, co, mo, be, b = build_materials(p)
    f3.temperature = 300.0
    g_f3, _, _ = build_geometry(p, f3, cl, co, mo, be, b)
    mod_f3 = openmc.Model(materials=mats_f3, geometry=g_f3, settings=build_settings(p))
    mod_f3.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_f3 = float(sp.keff.nominal_value)
        s_f3 = float(sp.keff.std_dev)
    print(f"  k(fuel 300 K)  = {k_f3:.5f} +/- {s_f3:.5f}")

    openmc.reset_auto_ids()
    mats_f15, f15, cl, co, mo, be, b = build_materials(p)
    f15.temperature = 1500.0
    g_f15, _, _ = build_geometry(p, f15, cl, co, mo, be, b)
    mod_f15 = openmc.Model(materials=mats_f15, geometry=g_f15, settings=build_settings(p))
    mod_f15.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_f15 = float(sp.keff.nominal_value)
        s_f15 = float(sp.keff.std_dev)
    print(f"  k(fuel 1500 K) = {k_f15:.5f} +/- {s_f15:.5f}")

    alpha_fuel = (((k_f15 - 1.0)/k_f15 - (k_f3 - 1.0)/k_f3) * 1e5) / 1200.0
    u_alpha_fuel = np.sqrt((s_f3/k_f3**2)**2 + (s_f15/k_f15**2)**2) / 1200.0 * 1e5
    print(f"  alpha_fuel = {alpha_fuel:+.2f} +/- {u_alpha_fuel:.2f} pcm/K (1-sigma)")

    # 2. Moderator Feedback: 300 K vs 900 K
    print("\n--- 2. Evaluating Moderator Feedback (300 K vs 900 K) ---")
    openmc.reset_auto_ids()
    mats_m3, f, cl, co, m3, be, b = build_materials(p)
    m3.temperature = 300.0
    g_m3, _, _ = build_geometry(p, f, cl, co, m3, be, b)
    mod_m3 = openmc.Model(materials=mats_m3, geometry=g_m3, settings=build_settings(p))
    mod_m3.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_m3 = float(sp.keff.nominal_value)
        s_m3 = float(sp.keff.std_dev)
    print(f"  k(mod 300 K) = {k_m3:.5f} +/- {s_m3:.5f}")

    openmc.reset_auto_ids()
    mats_m9, f, cl, co, m9, be, b = build_materials(p)
    m9.temperature = 900.0
    g_m9, _, _ = build_geometry(p, f, cl, co, m9, be, b)
    mod_m9 = openmc.Model(materials=mats_m9, geometry=g_m9, settings=build_settings(p))
    mod_m9.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_m9 = float(sp.keff.nominal_value)
        s_m9 = float(sp.keff.std_dev)
    print(f"  k(mod 900 K) = {k_m9:.5f} +/- {s_m9:.5f}")

    alpha_mod = (((k_m9 - 1.0)/k_m9 - (k_m3 - 1.0)/k_m3) * 1e5) / 600.0
    u_alpha_mod = np.sqrt((s_m3/k_m3**2)**2 + (s_m9/k_m9**2)**2) / 600.0 * 1e5
    print(f"  alpha_mod  = {alpha_mod:+.2f} +/- {u_alpha_mod:.2f} pcm/K (1-sigma)")

    # 3. Coolant Feedback: cold (900 K) vs hot (1300 K)
    print("\n--- 3. Evaluating Coolant Feedback (900 K vs 1300 K) ---")
    openmc.reset_auto_ids()
    p_cold = ReactorParameters(
        core_radius=p.core_radius, active_core_height=p.active_core_height,
        fuel_radius=p.fuel_radius, cladding_thickness=p.cladding_thickness,
        coolant_channel_thickness=p.coolant_channel_thickness, moderator_fuel_pitch=p.moderator_fuel_pitch,
        radial_reflector_thickness=p.radial_reflector_thickness, axial_reflector_thickness=p.axial_reflector_thickness,
        control_drum_radius=p.control_drum_radius, absorber_thickness=p.absorber_thickness,
        absorber_sector_angular_width=p.absorber_sector_angular_width, b10_enrichment=p.b10_enrichment,
        coolant_state='cold', particles=p.particles, batches=p.batches, inactive=p.inactive
    )
    mats_cc, f_c, cl_c, co_c, mo_c, be_c, b_c = build_materials(p_cold)
    g_cc, _, _ = build_geometry(p_cold, f_c, cl_c, co_c, mo_c, be_c, b_c)
    mod_cc = openmc.Model(materials=mats_cc, geometry=g_cc, settings=build_settings(p_cold))
    mod_cc.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_cc = float(sp.keff.nominal_value)
        s_cc = float(sp.keff.std_dev)
    print(f"  k(cool cold 900 K)  = {k_cc:.5f} +/- {s_cc:.5f}")

    openmc.reset_auto_ids()
    p_hot = ReactorParameters(
        core_radius=p.core_radius, active_core_height=p.active_core_height,
        fuel_radius=p.fuel_radius, cladding_thickness=p.cladding_thickness,
        coolant_channel_thickness=p.coolant_channel_thickness, moderator_fuel_pitch=p.moderator_fuel_pitch,
        radial_reflector_thickness=p.radial_reflector_thickness, axial_reflector_thickness=p.axial_reflector_thickness,
        control_drum_radius=p.control_drum_radius, absorber_thickness=p.absorber_thickness,
        absorber_sector_angular_width=p.absorber_sector_angular_width, b10_enrichment=p.b10_enrichment,
        coolant_state='hot', particles=p.particles, batches=p.batches, inactive=p.inactive
    )
    mats_ch, f_h, cl_h, co_h, mo_h, be_h, b_h = build_materials(p_hot)
    g_ch, _, _ = build_geometry(p_hot, f_h, cl_h, co_h, mo_h, be_h, b_h)
    mod_ch = openmc.Model(materials=mats_ch, geometry=g_ch, settings=build_settings(p_hot))
    mod_ch.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_ch = float(sp.keff.nominal_value)
        s_ch = float(sp.keff.std_dev)
    print(f"  k(cool hot 1300 K) = {k_ch:.5f} +/- {s_ch:.5f}")

    alpha_cool = (((k_ch - 1.0)/k_ch - (k_cc - 1.0)/k_cc) * 1e5) / 400.0
    u_alpha_cool = np.sqrt((s_cc/k_cc**2)**2 + (s_ch/k_ch**2)**2) / 400.0 * 1e5
    print(f"  alpha_cool = {alpha_cool:+.2f} +/- {u_alpha_cool:.2f} pcm/K (1-sigma)")

    alpha_tot = alpha_fuel + alpha_mod + alpha_cool
    u_tot = np.sqrt(u_alpha_fuel**2 + u_alpha_mod**2 + u_alpha_cool**2)

    print("\n" + "=" * 80)
    print("TOTAL TEMPERATURE COEFFICIENT SUMMARY:")
    print("=" * 80)
    print(f"  alpha_fuel : {alpha_fuel:+.2f} +/- {u_alpha_fuel:.2f} pcm/K (1-sigma), +/- {2*u_alpha_fuel:.2f} (2-sigma), +/- {3*u_alpha_fuel:.2f} (3-sigma)")
    print(f"  alpha_mod  : {alpha_mod:+.2f} +/- {u_alpha_mod:.2f} pcm/K (1-sigma), +/- {2*u_alpha_mod:.2f} (2-sigma), +/- {3*u_alpha_mod:.2f} (3-sigma)")
    print(f"  alpha_cool : {alpha_cool:+.2f} +/- {u_alpha_cool:.2f} pcm/K (1-sigma), +/- {2*u_alpha_cool:.2f} (2-sigma), +/- {3*u_alpha_cool:.2f} (3-sigma)")
    print(f"  alpha_total: {alpha_tot:+.2f} +/- {u_tot:.2f} pcm/K (1-sigma), +/- {2*u_tot:.2f} (2-sigma), +/- {3*u_tot:.2f} (3-sigma)")
    print(f"  3-sigma upper bound for total TMC: {alpha_tot + 3*u_tot:+.2f} pcm/K")

    os.chdir(orig_dir)
    return alpha_fuel, u_alpha_fuel, alpha_mod, u_alpha_mod, alpha_cool, u_alpha_cool, alpha_tot, u_tot

if __name__ == "__main__":
    run_tmc_test()
