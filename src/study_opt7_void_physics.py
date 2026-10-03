"""
Comprehensive Void Physics & Moderation Analysis for Opt-7 Candidate Development
================================================================================
Investigates the microscopic and macroscopic physics governing void reactivity:
1. Verification of the 'Under-Moderated' status:
   - Evaluates k_inf vs pitch (1.45 cm to 1.80 cm) to determine peak of k_inf curve
   - Confirms sign of d(k_inf)/d(V_mod)
2. Energy spectrum & reaction rate breakdown:
   - Thermal (< 0.625 eV), Epithermal (0.625 eV - 100 keV), Fast (> 100 keV)
   - Quantifies Xenon-129 and Xenon-131 absorption rates vs U-235 fission
3. Coolant composition sweep with high statistics:
   - He-Xe 72/28 (M = 40.0 g/mol)
   - He-Xe 80/20 (M = 29.5 g/mol)
   - He-Xe 85/15 (M = 23.1 g/mol)
   - He-Xe 90/10 (M = 16.7 g/mol)
   - He-Xe 95/5  (M = 10.4 g/mol)
   - Pure Helium (M = 4.0 g/mol)
4. Coolant channel gap sweep:
   - t_cool = 0.045 cm, 0.040 cm, 0.035 cm
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

def run_kinf_moderation_curve(work_dir="scratch/study_opt7_kinf"):
    """Evaluates infinite lattice k_inf vs pitch to establish the under-moderated regime."""
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    print("=" * 80)
    print("PART 1: INFINITE PIN-CELL LATTICE k_inf vs PITCH (MODERATION CURVE)")
    print("=" * 80)

    pitches = [1.45, 1.50, 1.55, 1.615, 1.68, 1.75, 1.85, 2.00]
    results = []

    for p_val in pitches:
        openmc.reset_auto_ids()
        # Single pin cell with reflective boundary
        rp = ReactorParameters(
            fuel_radius=0.58, cladding_thickness=0.05, coolant_channel_thickness=0.045,
            moderator_fuel_pitch=p_val, particles=8000, batches=35, inactive=10
        )
        mats, fuel, clad, cool, mod, beo, b4c = build_materials(rp)

        s_fuel = openmc.ZCylinder(r=rp.fuel_radius)
        s_clad = openmc.ZCylinder(r=rp.clad_outer_radius)
        s_chan = openmc.ZCylinder(r=rp.channel_outer_radius)

        # Reflective box for unit cell
        hp = p_val / 2.0
        x_min = openmc.XPlane(x0=-hp, boundary_type='reflective')
        x_max = openmc.XPlane(x0=+hp, boundary_type='reflective')
        y_min = openmc.YPlane(y0=-hp, boundary_type='reflective')
        y_max = openmc.YPlane(y0=+hp, boundary_type='reflective')
        z_min = openmc.ZPlane(z0=-10.0, boundary_type='reflective')
        z_max = openmc.ZPlane(z0=+10.0, boundary_type='reflective')

        c_fuel = openmc.Cell(fill=fuel, region=-s_fuel & +x_min & -x_max & +y_min & -y_max & +z_min & -z_max)
        c_clad = openmc.Cell(fill=clad, region=+s_fuel & -s_clad & +x_min & -x_max & +y_min & -y_max & +z_min & -z_max)
        c_cool = openmc.Cell(fill=cool, region=+s_clad & -s_chan & +x_min & -x_max & +y_min & -y_max & +z_min & -z_max)
        c_mod = openmc.Cell(fill=mod, region=+s_chan & +x_min & -x_max & +y_min & -y_max & +z_min & -z_max)

        geom = openmc.Geometry([c_fuel, c_clad, c_cool, c_mod])
        sett = openmc.Settings()
        sett.run_mode = 'eigenvalue'
        sett.particles = rp.particles
        sett.batches = rp.batches
        sett.inactive = rp.inactive
        sett.source = openmc.IndependentSource(space=openmc.stats.Box([-hp*0.5, -hp*0.5, -5], [hp*0.5, hp*0.5, 5]))
        sett.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}

        mod_pin = openmc.Model(materials=mats, geometry=geom, settings=sett)
        mod_pin.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{rp.batches}.h5") as sp:
            k_inf = float(sp.keff.nominal_value)
            s_inf = float(sp.keff.std_dev)

        v_f = np.pi * rp.fuel_radius**2
        v_mod = p_val**2 - np.pi * rp.channel_outer_radius**2
        ratio = v_mod / v_f

        rec = {"pitch": p_val, "v_mod_ratio": ratio, "k_inf": k_inf, "s_inf": s_inf}
        results.append(rec)
        print(f"Pitch = {p_val:5.3f} cm | V_mod/V_fuel = {ratio:5.3f} | k_inf = {k_inf:.5f} +/- {s_inf:.5f}")

    os.chdir(orig_dir)
    csv_path = "results/opt7/kinf_moderation_curve.csv"
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Pitch_cm", "Vmod_Vfuel_cell", "kinf", "std_dev"])
        for r in results:
            writer.writerow([f"{r['pitch']:.3f}", f"{r['v_mod_ratio']:.3f}", f"{r['k_inf']:.5f}", f"{r['s_inf']:.5f}"])
    print(f"Saved: {csv_path}")
    return results

def run_coolant_composition_and_void_matrix(work_dir="scratch/study_opt7_void"):
    """Evaluates void reactivity and reaction rate decomposition across gas mixtures."""
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    print("\n" + "=" * 80)
    print("PART 2: SYSTEMATIC COOLANT COMPOSITION & VOID REACTIVITY (HIGH STATS)")
    print("=" * 80)

    # Base Opt-6 geometry
    p = ReactorParameters(
        core_radius=17.90, active_core_height=35.00,
        fuel_radius=0.58, cladding_thickness=0.05,
        coolant_channel_thickness=0.045, moderator_fuel_pitch=1.615,
        radial_reflector_thickness=9.60, axial_reflector_thickness=7.50,
        control_drum_radius=4.30, absorber_thickness=1.45,
        absorber_sector_angular_width=140.0, b10_enrichment=0.90,
        particles=20000, batches=50, inactive=10
    )

    compositions = [
        {"name": "He-Xe 72/28 (M=40.0)", "he": 0.72, "xe": 0.28, "M": 40.038},
        {"name": "He-Xe 80/20 (M=29.5)", "he": 0.80, "xe": 0.20, "M": 29.461},
        {"name": "He-Xe 85/15 (M=23.1)", "he": 0.85, "xe": 0.15, "M": 23.096},
        {"name": "He-Xe 90/10 (M=16.7)", "he": 0.90, "xe": 0.10, "M": 16.732},
        {"name": "He-Xe 95/5  (M=10.4)", "he": 0.95, "xe": 0.05, "M": 10.367},
        {"name": "Pure He    (M=4.0)",  "he": 1.00, "xe": 0.00, "M":  4.003},
    ]

    # 1. Voided State (Run with very high statistics: 25k particles, 50 batches)
    print("\n--- 1. Evaluating Voided Reference Core (Vacuum) ---")
    openmc.reset_auto_ids()
    energy_filter_v = openmc.EnergyFilter([0.0, 0.625, 100.0e3, 20.0e6])
    p_vd = ReactorParameters(
        core_radius=p.core_radius, active_core_height=p.active_core_height,
        fuel_radius=p.fuel_radius, cladding_thickness=p.cladding_thickness,
        coolant_channel_thickness=p.coolant_channel_thickness,
        moderator_fuel_pitch=p.moderator_fuel_pitch,
        radial_reflector_thickness=p.radial_reflector_thickness,
        axial_reflector_thickness=p.axial_reflector_thickness,
        control_drum_radius=p.control_drum_radius,
        absorber_thickness=p.absorber_thickness,
        absorber_sector_angular_width=p.absorber_sector_angular_width,
        b10_enrichment=p.b10_enrichment,
        coolant_state='voided',
        particles=25000, batches=50, inactive=10
    )
    p_vd.coolant_density = 1.0e-7
    mats_v, f_v, cl_v, co_v, mo_v, be_v, b_v = build_materials(p_vd)
    g_v, _, _ = build_geometry(p_vd, f_v, cl_v, co_v, mo_v, be_v, b_v)
    t_v = openmc.Tallies()
    t_spec_v = openmc.Tally(name="spectrum_void")
    t_spec_v.filters = [energy_filter_v]
    t_spec_v.scores = ["flux", "fission", "absorption"]
    t_v.append(t_spec_v)
    s_v = build_settings(p_vd)
    mod_v = openmc.Model(materials=mats_v, geometry=g_v, settings=s_v, tallies=t_v)
    mod_v.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p_vd.batches}.h5") as sp:
        k_void = float(sp.keff.nominal_value)
        s_void = float(sp.keff.std_dev)
    print(f"Voided Reference: k_void = {k_void:.5f} +/- {s_void:.5f}")

    comp_results = []
    for c in compositions:
        openmc.reset_auto_ids()
        mats, fuel, clad, cool, mod, beo, b4c = build_materials(p)
        # Custom coolant definition
        cool.nuclides.clear()
        # Ideal gas law at 2.0 MPa, 1050 K: rho = 0.0002291 * M [g/cm3]
        rho_calc = 0.0002291 * c["M"]
        cool.set_density("g/cm3", rho_calc)
        if c["he"] > 0:
            cool.add_nuclide("He4", c["he"], "ao")
        if c["xe"] > 0:
            cool.add_element("Xe", c["xe"], "ao")

        g_op, _, _ = build_geometry(p, fuel, clad, cool, mod, beo, b4c)

        # Tallies to get flux spectrum and absorption in coolant
        energy_filter = openmc.EnergyFilter([0.0, 0.625, 100.0e3, 20.0e6])
        tallies = openmc.Tallies()
        t_spec = openmc.Tally(name="spectrum_operating")
        t_spec.filters = [energy_filter]
        t_spec.scores = ["flux", "fission", "absorption"]
        tallies.append(t_spec)

        # Coolant absorption tally
        cool_filter = openmc.MaterialFilter(cool)
        t_cool_abs = openmc.Tally(name="coolant_absorption")
        t_cool_abs.filters = [cool_filter]
        t_cool_abs.scores = ["absorption"]
        tallies.append(t_cool_abs)

        settings = build_settings(p)
        model = openmc.Model(materials=mats, geometry=g_op, settings=settings, tallies=tallies)
        model.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
            k_op = float(sp.keff.nominal_value)
            s_op = float(sp.keff.std_dev)
            t_sp = sp.get_tally(name="spectrum_operating")
            flx = t_sp.get_values(scores=["flux"]).flatten()
            tot_flx = np.sum(flx)
            f_therm = flx[0] / tot_flx if tot_flx > 0 else 0
            f_epi   = flx[1] / tot_flx if tot_flx > 0 else 0
            f_fast  = flx[2] / tot_flx if tot_flx > 0 else 0

            t_ca = sp.get_tally(name="coolant_absorption")
            cool_abs_val = float(t_ca.get_values(scores=["absorption"]).flatten()[0])

        # Exact reactivity formulation:
        # Delta-rho = (1/k_nom - 1/k_void) * 1e5
        drho_void = (1.0 / k_op - 1.0 / k_void) * 1e5
        u_void = np.sqrt((s_op/k_op**2)**2 + (s_void/k_void**2)**2) * 1e5
        bound_3sigma = drho_void + 3.0 * u_void

        rec = {
            "name": c["name"], "he": c["he"], "xe": c["xe"], "M": c["M"], "density": rho_calc,
            "k_op": k_op, "s_op": s_op, "k_void": k_void, "s_void": s_void,
            "drho_void": drho_void, "u_void": u_void, "bound_3sigma": bound_3sigma,
            "f_therm": f_therm, "f_epi": f_epi, "f_fast": f_fast,
            "cool_abs": cool_abs_val
        }
        comp_results.append(rec)
        print(f"\n{c['name']:24s} | rho = {rho_calc:.5f} g/cm3")
        print(f"  k_op = {k_op:.5f} +/- {s_op:.5f} | k_void = {k_void:.5f} +/- {s_void:.5f}")
        print(f"  Delta-rho_void = {drho_void:+.1f} +/- {u_void:.1f} pcm (3-sigma bound: {bound_3sigma:+.1f} pcm)")
        print(f"  Spectrum: Therm = {f_therm*100:.1f}%, Epi = {f_epi*100:.1f}%, Fast = {f_fast*100:.1f}% | Coolant Abs = {cool_abs_val:.6f}")

    os.chdir(orig_dir)
    csv_path = "results/opt7/coolant_composition_and_void_physics.csv"
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Composition", "He_ao", "Xe_ao", "M_g_mol", "Density_g_cm3",
            "k_op", "s_op", "k_void", "s_void",
            "Delta_rho_void_pcm", "u_void_pcm", "Bound_3sigma_pcm",
            "Flux_Frac_Thermal", "Flux_Frac_Epithermal", "Flux_Frac_Fast", "Coolant_Absorption_Rate"
        ])
        for r in comp_results:
            writer.writerow([
                r["name"], f"{r['he']:.2f}", f"{r['xe']:.2f}", f"{r['M']:.3f}", f"{r['density']:.5f}",
                f"{r['k_op']:.5f}", f"{r['s_op']:.5f}", f"{r['k_void']:.5f}", f"{r['s_void']:.5f}",
                f"{r['drho_void']:+.1f}", f"{r['u_void']:.1f}", f"{r['bound_3sigma']:+.1f}",
                f"{r['f_therm']:.4f}", f"{r['f_epi']:.4f}", f"{r['f_fast']:.4f}", f"{r['cool_abs']:.6e}"
            ])
    print(f"\nSaved physics dataset to: {csv_path}")
    return comp_results

if __name__ == "__main__":
    if not os.path.exists("results/opt7/kinf_moderation_curve.csv"):
        run_kinf_moderation_curve()
    run_coolant_composition_and_void_matrix()
