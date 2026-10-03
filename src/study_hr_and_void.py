"""
Systematic Study of Core H/R Ratio, Moderator Fraction, and Coolant Physics
==========================================================================
Investigates the primary physics levers identified for Opt-6:
1. Core H/R ratio sweep at constant core volume (V = 35.0 L):
   H/R = 1.4, 1.6, 1.85, 2.0, 2.2, 2.5
2. Coolant void reactivity decomposition:
   - Evaluates Xe absorption removal vs He moderation/leakage
   - Tests coolant mixtures: He-Xe 72/28 (M=40), He-Xe 85/15 (M=23), He-Xe 90/10 (M=17), pure He (M=4)
3. Moderator ratio / under-moderation sweep (Ade et al. 2022):
   - Pitch: 1.635 cm, 1.55 cm, 1.48 cm
   - Fuel radius: 0.58 cm vs 0.62 cm
   - Evaluates MTC, Doppler, Total TMC, and Void reactivity
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

def run_hr_sensitivity(work_dir="scratch/study_opt6"):
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    print("=" * 80)
    print("PART 1: SYSTEMATIC H/R GEOMETRIC SENSITIVITY STUDY (CONSTANT V = 35.0 L)")
    print("=" * 80)

    # Constant core volume: V = pi * R^2 * H = 35,000 cm3
    # H = (H/R) * R -> pi * R^3 * (H/R) = V -> R = (V / (pi * (H/R)))^(1/3)
    V_core = 35000.0  # cm3
    hr_targets = [1.4, 1.6, 1.85, 2.0, 2.2, 2.5, 2.72]

    results = []
    for hr in hr_targets:
        R_c = (V_core / (np.pi * hr))**(1.0 / 3.0)
        H_c = hr * R_c
        hd = hr / 2.0  # H/D

        p = ReactorParameters(
            core_radius=R_c,
            active_core_height=H_c,
            moderator_fuel_pitch=1.635,
            radial_reflector_thickness=9.60,
            axial_reflector_thickness=7.60,
            control_drum_radius=4.30,
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

        openmc.reset_auto_ids()
        p.drum_angle_deg = 0.0
        mats, f, cl, co, m, be, b = build_materials(p)
        geom, _, _ = build_geometry(p, f, cl, co, m, be, b)
        sett = build_settings(p)
        mod = openmc.Model(materials=mats, geometry=geom, settings=sett)
        mod.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
            k_out = float(sp.keff.nominal_value)
            s_out = float(sp.keff.std_dev)

        rec = {
            "hr": hr, "hd": hd, "r_core": R_c, "h_core": H_c, "pins": pins,
            "fuel_mass_kg": inv["total_fuel_mass_kg"],
            "core_mass_kg": inv['items'][0]['mass_kg'] + inv['items'][1]['mass_kg'] + inv['items'][2]['mass_kg'] + inv['items'][3]['mass_kg'],
            "refl_mass_kg": inv['total_reflector_mass_kg'] + inv['total_absorber_mass_kg'],
            "total_mass_kg": inv["total_mass_kg"],
            "total_volume_L": inv["total_volume_L"],
            "k_out": k_out, "s_out": s_out
        }
        results.append(rec)
        print(f"H/R = {hr:4.2f} (H/D = {hd:4.2f}) | R = {R_c:5.2f} cm, H = {H_c:5.2f} cm | Pins: {pins:3d} | Mass: {inv['total_mass_kg']:6.2f} kg | k_out = {k_out:.5f} +/- {s_out:.5f}")

    os.chdir(orig_dir)
    csv_path = "results/opt6/hr_sensitivity_study.csv"
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["H/R", "H/D", "R_core [cm]", "H_core [cm]", "Pins", "Fuel Mass [kg]", "Core Mass [kg]", "Refl Mass [kg]", "Total Mass [kg]", "Total Volume [L]", "k_out", "s_out"])
        for r in results:
            writer.writerow([
                f"{r['hr']:.2f}", f"{r['hd']:.2f}", f"{r['r_core']:.2f}", f"{r['h_core']:.2f}",
                r["pins"], f"{r['fuel_mass_kg']:.2f}", f"{r['core_mass_kg']:.2f}", f"{r['refl_mass_kg']:.2f}",
                f"{r['total_mass_kg']:.2f}", f"{r['total_volume_L']:.2f}", f"{r['k_out']:.5f}", f"{r['s_out']:.5f}"
            ])
    print(f"Saved: {csv_path}")
    return results

def run_moderation_and_void_study(work_dir="scratch/study_opt6"):
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    print("\n" + "=" * 80)
    print("PART 2: MODERATION & COOLANT VOID REACTIVITY DECOMPOSITION STUDY")
    print("=" * 80)

    # Test cases:
    # 1. Base Opt-5 (P = 1.635 cm, r_fuel = 0.58 cm, V_mod/V_fuel = 1.15)
    # 2. Under-moderated A (P = 1.55 cm, r_fuel = 0.58 cm, V_mod/V_fuel = 0.95)
    # 3. Under-moderated B (P = 1.50 cm, r_fuel = 0.60 cm, V_mod/V_fuel = 0.73)
    # 4. Under-moderated C (P = 1.46 cm, r_fuel = 0.62 cm, V_mod/V_fuel = 0.58)

    cases = [
        {"name": "Opt-5 Base (P=1.635, rf=0.58)", "pitch": 1.635, "rf": 0.58, "gap": 0.05, "rc": 16.0, "hc": 43.5},
        {"name": "Under-Mod A (P=1.55, rf=0.58)", "pitch": 1.550, "rf": 0.58, "gap": 0.05, "rc": 15.5, "hc": 43.5},
        {"name": "Under-Mod B (P=1.52, rf=0.60)", "pitch": 1.520, "rf": 0.60, "gap": 0.045, "rc": 15.2, "hc": 43.5},
        {"name": "Under-Mod C (P=1.48, rf=0.62)", "pitch": 1.480, "rf": 0.62, "gap": 0.040, "rc": 15.0, "hc": 43.5},
    ]

    coolant_types = [
        {"type": "operating", "desc": "He-Xe 72/28 (M=40, 2.0 MPa)"},
        {"type": "voided", "desc": "Voided (0 MPa vacuum)"}
    ]

    mod_results = []
    for c in cases:
        p_base = ReactorParameters(
            fuel_radius=c["rf"],
            cladding_thickness=0.05,
            coolant_channel_thickness=c["gap"],
            moderator_fuel_pitch=c["pitch"],
            core_radius=c["rc"],
            active_core_height=c["hc"],
            radial_reflector_thickness=9.60,
            axial_reflector_thickness=7.60,
            control_drum_radius=4.30,
            absorber_thickness=1.45,
            absorber_sector_angular_width=140.0,
            particles=6000,
            batches=32,
            inactive=8
        )
        pins = count_pins(p_base)
        inv = calculate_inventory(p_base, pins)

        # 1. Operating keff
        openmc.reset_auto_ids()
        p_base.coolant_state = 'operating'
        mats, f, cl, co, m, be, b = build_materials(p_base)
        geom, _, _ = build_geometry(p_base, f, cl, co, m, be, b)
        sett = build_settings(p_base)
        mod_op = openmc.Model(materials=mats, geometry=geom, settings=sett)
        mod_op.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_base.batches}.h5") as sp:
            k_op = float(sp.keff.nominal_value)
            s_op = float(sp.keff.std_dev)

        # 2. Voided keff
        openmc.reset_auto_ids()
        p_base.coolant_state = 'voided'
        mats_v, f_v, cl_v, co_v, m_v, be_v, b_v = build_materials(p_base)
        geom_v, _, _ = build_geometry(p_base, f_v, cl_v, co_v, m_v, be_v, b_v)
        sett_v = build_settings(p_base)
        mod_vd = openmc.Model(materials=mats_v, geometry=geom_v, settings=sett_v)
        mod_vd.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_base.batches}.h5") as sp:
            k_vd = float(sp.keff.nominal_value)
            s_vd = float(sp.keff.std_dev)

        # 3. Doppler and MTC check
        # Fuel at 1200 K
        openmc.reset_auto_ids()
        p_base.coolant_state = 'operating'
        mats_d, f_d, cl_d, co_d, m_d, be_d, b_d = build_materials(p_base)
        f_d.temperature = 1200.0
        geom_d, _, _ = build_geometry(p_base, f_d, cl_d, co_d, m_d, be_d, b_d)
        sett_d = build_settings(p_base)
        mod_dp = openmc.Model(materials=mats_d, geometry=geom_d, settings=sett_d)
        mod_dp.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_base.batches}.h5") as sp:
            k_dop = float(sp.keff.nominal_value)
            s_dop = float(sp.keff.std_dev)

        # Moderator at 900 K
        openmc.reset_auto_ids()
        mats_m, f_m, cl_m, co_m, m_m, be_m, b_m = build_materials(p_base)
        m_m.temperature = 900.0
        geom_m, _, _ = build_geometry(p_base, f_m, cl_m, co_m, m_m, be_m, b_m)
        sett_m = build_settings(p_base)
        mod_mt = openmc.Model(materials=mats_m, geometry=geom_m, settings=sett_m)
        mod_mt.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_base.batches}.h5") as sp:
            k_mtc = float(sp.keff.nominal_value)
            s_mtc = float(sp.keff.std_dev)

        drho_void = ((k_vd - 1.0)/k_vd - (k_op - 1.0)/k_op) * 1e5
        u_void = np.sqrt((s_vd/k_vd**2)**2 + (s_op/k_op**2)**2) * 1e5

        alpha_fuel = (((k_dop - 1.0)/k_dop - (k_op - 1.0)/k_op) * 1e5) / (1200.0 - 1050.0) if 1200.0 != 1050.0 else 0.0
        u_alpha_fuel = np.sqrt((s_dop/k_dop**2)**2 + (s_op/k_op**2)**2) / 150.0 * 1e5

        alpha_mod = (((k_mtc - 1.0)/k_mtc - (k_op - 1.0)/k_op) * 1e5) / (900.0 - 1050.0)
        u_alpha_mod = np.sqrt((s_mtc/k_mtc**2)**2 + (s_op/k_op**2)**2) / 150.0 * 1e5

        # V_mod / V_fuel calculation
        v_fuel_pin = np.pi * c["rf"]**2
        v_cell = c["pitch"]**2
        v_channel = np.pi * (c["rf"] + 0.05 + c["gap"])**2
        v_mod_cell = v_cell - v_channel
        v_ratio = v_mod_cell / v_fuel_pin

        r_entry = {
            "name": c["name"], "pitch": c["pitch"], "rf": c["rf"], "gap": c["gap"],
            "v_ratio": v_ratio, "pins": pins, "mass_kg": inv["total_mass_kg"],
            "k_op": k_op, "s_op": s_op,
            "k_vd": k_vd, "s_vd": s_vd,
            "drho_void": drho_void, "u_void": u_void,
            "alpha_fuel": alpha_fuel, "u_alpha_fuel": u_alpha_fuel,
            "alpha_mod": alpha_mod, "u_alpha_mod": u_alpha_mod
        }
        mod_results.append(r_entry)
        print(f"\nCase: {c['name']} (V_mod/V_fuel = {v_ratio:.2f})")
        print(f"  Pins: {pins} | Mass: {inv['total_mass_kg']:.2f} kg | k_op = {k_op:.5f}")
        print(f"  VOID REACTIVITY : {drho_void:+.1f} +/- {u_void:.1f} pcm")
        print(f"  FUEL DOPPLER    : {alpha_fuel:+.2f} +/- {u_alpha_fuel:.2f} pcm/K")
        print(f"  MODERATOR COEFF : {alpha_mod:+.2f} +/- {u_alpha_mod:.2f} pcm/K")

    os.chdir(orig_dir)
    csv_path = "results/opt6/moderation_and_void_study.csv"
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Case", "Pitch [cm]", "r_fuel [cm]", "Gap [cm]", "Vmod/Vfuel", "Pins", "Mass [kg]", "k_op", "s_op", "k_void", "s_void", "Void drho [pcm]", "u_void [pcm]", "alpha_fuel [pcm/K]", "u_alpha_fuel [pcm/K]", "alpha_mod [pcm/K]", "u_alpha_mod [pcm/K]"])
        for r in mod_results:
            writer.writerow([
                r["name"], f"{r['pitch']:.3f}", f"{r['rf']:.2f}", f"{r['gap']:.3f}",
                f"{r['v_ratio']:.2f}", r["pins"], f"{r['mass_kg']:.2f}",
                f"{r['k_op']:.5f}", f"{r['s_op']:.5f}", f"{r['k_vd']:.5f}", f"{r['s_vd']:.5f}",
                f"{r['drho_void']:+.1f}", f"{r['u_void']:.1f}",
                f"{r['alpha_fuel']:+.2f}", f"{r['u_alpha_fuel']:.2f}",
                f"{r['alpha_mod']:+.2f}", f"{r['u_alpha_mod']:.2f}"
            ])
    print(f"Saved: {csv_path}")
    return mod_results

if __name__ == "__main__":
    run_hr_sensitivity()
    run_moderation_and_void_study()
