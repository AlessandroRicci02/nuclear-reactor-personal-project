"""
Parametric Physics Exploration for Robust Negative Void Reactivity in Opt-7
===========================================================================
Investigates the interplay between:
1. He-Xe composition: 72/28, 85/15, 90/10, 95/5, 98/2, Pure He
2. Coolant channel gap: 0.035 cm, 0.045 cm, 0.055 cm
3. Operating pressure: 1.2 MPa, 1.6 MPa, 2.0 MPa
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

def run_simulation(p, rho_cool, he_frac, xe_frac):
    openmc.reset_auto_ids()
    mats, f, cl, co, mo, be, b = build_materials(p)
    co.nuclides.clear()
    co.set_density("g/cm3", rho_cool)
    co.add_nuclide("He4", he_frac, "ao")
    if xe_frac > 0:
        co.add_element("Xe", xe_frac, "ao")
    geom, pins, _ = build_geometry(p, f, cl, co, mo, be, b)
    sett = build_settings(p)
    mod = openmc.Model(materials=mats, geometry=geom, settings=sett)
    mod.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_val = float(sp.keff.nominal_value)
        k_std = float(sp.keff.std_dev)
    return k_val, k_std, len(pins)

def main():
    work_dir = "scratch/parametric_void_study"
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    print("=" * 80)
    print("PARAMETRIC VOID REACTIVITY OPTIMIZATION & SENSITIVITY FOR OPT-7")
    print("=" * 80)

    # Base Candidate Opt7-D geometry: R=17.80 cm, H=34.80 cm, Pitch=1.612 cm
    p_base = ReactorParameters(
        core_radius=17.80, active_core_height=34.80,
        fuel_radius=0.58, cladding_thickness=0.05,
        moderator_fuel_pitch=1.612,
        radial_reflector_thickness=9.60, axial_reflector_thickness=7.50,
        control_drum_radius=4.30, absorber_thickness=1.45,
        absorber_sector_angular_width=140.0, b10_enrichment=0.90,
        particles=40000, batches=60, inactive=12
    )

    # Dictionary of voided reference states by gap
    void_refs = {}

    gaps = [0.045, 0.035, 0.055]
    for g in gaps:
        print(f"\n--- Simulating Voided Reference State for Gap = {g*10:.2f} mm ---")
        p_v = ReactorParameters(
            core_radius=p_base.core_radius, active_core_height=p_base.active_core_height,
            fuel_radius=p_base.fuel_radius, cladding_thickness=p_base.cladding_thickness,
            coolant_channel_thickness=g, moderator_fuel_pitch=p_base.moderator_fuel_pitch,
            radial_reflector_thickness=p_base.radial_reflector_thickness,
            axial_reflector_thickness=p_base.axial_reflector_thickness,
            control_drum_radius=p_base.control_drum_radius,
            absorber_thickness=p_base.absorber_thickness,
            absorber_sector_angular_width=p_base.absorber_sector_angular_width,
            b10_enrichment=p_base.b10_enrichment,
            particles=60000, batches=70, inactive=15  # High statistics reference
        )
        kv, sv, n_pins = run_simulation(p_v, 1.0e-7, 1.0, 0.0)
        void_refs[g] = (kv, sv)
        print(f"  Voided Reference (gap {g*10:.2f} mm): k_void = {kv:.5f} +/- {sv:.5f}")

    cases = [
        # 1. Systematic gas sweep at 2.0 MPa, 0.45 mm gap
        ("He-Xe 72/28, 2.0 MPa",  0.72, 0.28, 2.0, 0.045),
        ("He-Xe 80/20, 2.0 MPa",  0.80, 0.20, 2.0, 0.045),
        ("He-Xe 85/15, 2.0 MPa",  0.85, 0.15, 2.0, 0.045),
        ("He-Xe 90/10, 2.0 MPa",  0.90, 0.10, 2.0, 0.045),
        ("He-Xe 95/5,  2.0 MPa",  0.95, 0.05, 2.0, 0.045),
        ("He-Xe 98/2,  2.0 MPa",  0.98, 0.02, 2.0, 0.045),
        ("Pure He,     2.0 MPa",  1.00, 0.00, 2.0, 0.045),

        # 2. Pressure sweep for He-Xe 90/10 and 95/5
        ("He-Xe 90/10, 1.5 MPa",  0.90, 0.10, 1.5, 0.045),
        ("He-Xe 90/10, 1.2 MPa",  0.90, 0.10, 1.2, 0.045),
        ("He-Xe 95/5,  1.5 MPa",  0.95, 0.05, 1.5, 0.045),
        ("He-Xe 95/5,  1.2 MPa",  0.95, 0.05, 1.2, 0.045),

        # 3. Gap sweep for He-Xe 95/5 and 90/10
        ("He-Xe 95/5, gap 0.35mm", 0.95, 0.05, 2.0, 0.035),
        ("He-Xe 95/5, gap 0.55mm", 0.95, 0.05, 2.0, 0.055),
        ("He-Xe 90/10, gap 0.35mm",0.90, 0.10, 2.0, 0.035),
    ]

    results = []
    print("\n--- Evaluating Operating States ---")
    for label, he_frac, xe_frac, press, gap in cases:
        p_c = ReactorParameters(
            core_radius=p_base.core_radius, active_core_height=p_base.active_core_height,
            fuel_radius=p_base.fuel_radius, cladding_thickness=p_base.cladding_thickness,
            coolant_channel_thickness=gap, moderator_fuel_pitch=p_base.moderator_fuel_pitch,
            radial_reflector_thickness=p_base.radial_reflector_thickness,
            axial_reflector_thickness=p_base.axial_reflector_thickness,
            control_drum_radius=p_base.control_drum_radius,
            absorber_thickness=p_base.absorber_thickness,
            absorber_sector_angular_width=p_base.absorber_sector_angular_width,
            b10_enrichment=p_base.b10_enrichment,
            particles=p_base.particles, batches=p_base.batches, inactive=p_base.inactive
        )
        M = he_frac * 4.002602 + xe_frac * 131.293
        rho_op = (press * 1.0e6 * M * 1.0e-3) / (8.31446 * 1050.0) * 1.0e-3

        k_op, s_op, pins_count = run_simulation(p_c, rho_op, he_frac, xe_frac)
        k_vd, s_vd = void_refs[gap]

        drho = (1.0 / k_op - 1.0 / k_vd) * 1.0e5
        u_drho = np.sqrt((s_op / k_op**2)**2 + (s_vd / k_vd**2)**2) * 1.0e5
        bound_3sigma = drho + 3.0 * u_drho

        inv = calculate_inventory(p_c, pins_count)
        mass_kg = inv["total_mass_kg"]
        vol_L = inv["total_volume_L"]

        rec = {
            "label": label, "he": he_frac, "xe": xe_frac, "M": M, "pressure_MPa": press,
            "gap_cm": gap, "k_op": k_op, "s_op": s_op, "k_vd": k_vd, "s_vd": s_vd,
            "drho": drho, "u_drho": u_drho, "bound_3sigma": bound_3sigma,
            "mass_kg": mass_kg, "vol_L": vol_L, "pins": pins_count
        }
        results.append(rec)
        print(f"[{label:24s}] k_op={k_op:.5f}+/-{s_op:.5f} | k_vd={k_vd:.5f}+/-{s_vd:.5f} | drho={drho:+6.1f}+/-{u_drho:4.1f} pcm | 3-sig bound={bound_3sigma:+6.1f} pcm | M={mass_kg:.1f}kg")

    os.chdir(orig_dir)

    out_csv = "results/opt7/void/parametric_void_study.csv"
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Configuration", "He_ao", "Xe_ao", "M_g_mol", "Pressure_MPa", "Gap_cm", "Pins",
                         "Mass_kg", "Volume_L", "k_op", "s_op", "k_vd", "s_vd",
                         "Delta_rho_void_pcm", "u_drho_pcm", "Bound_3sigma_pcm"])
        for r in results:
            writer.writerow([
                r["label"], f"{r['he']:.2f}", f"{r['xe']:.2f}", f"{r['M']:.3f}",
                f"{r['pressure_MPa']:.2f}", f"{r['gap_cm']:.4f}", r["pins"],
                f"{r['mass_kg']:.2f}", f"{r['vol_L']:.2f}",
                f"{r['k_op']:.5f}", f"{r['s_op']:.5f}", f"{r['k_vd']:.5f}", f"{r['s_vd']:.5f}",
                f"{r['drho']:+.1f}", f"{r['u_drho']:.1f}", f"{r['bound_3sigma']:+.1f}"
            ])
    print(f"\nSaved comprehensive parametric study to {out_csv}")

if __name__ == "__main__":
    main()
