"""
Systematic Study of Coolant Composition and Void Reactivity Physics
===================================================================
Investigates the mechanisms governing the coolant void reactivity in NEP cores:
1. Coolant composition variations:
   - He-Xe 72/28 (M = 40.0 g/mol, standard heavy Brayton gas)
   - He-Xe 80/20 (M = 29.5 g/mol)
   - He-Xe 85/15 (M = 23.1 g/mol)
   - He-Xe 90/10 (M = 16.7 g/mol)
   - Pure He (M = 4.0 g/mol)
2. Spectrum hardening / under-moderation lever:
   - Tests how shifting V_mod / V_fuel from 1.15 down to 0.75-0.90 affects:
     * Coolant void reactivity (pcm)
     * Moderator temperature coefficient (pcm/K)
     * Fuel Doppler coefficient (pcm/K)
3. H/R ratio integration:
   - Evaluates whether optimal H/R (1.85) enhances axial leakage upon voiding.
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

def build_custom_materials(p: ReactorParameters, he_frac=0.72, xe_frac=0.28, coolant_density=None,
                           fuel_temp=None, mod_temp=None):
    """Build materials with custom coolant gas composition, density, and temperatures."""
    materials, fuel, cladding, coolant, moderator, reflector, absorber = build_materials(p)

    if fuel_temp is not None:
        fuel.temperature = fuel_temp
    if mod_temp is not None:
        moderator.temperature = mod_temp

    # Re-define coolant with custom composition
    coolant.nuclides.clear()
    m_mol = he_frac * 4.0026 + xe_frac * 131.293
    # Ideal gas law at 2.0 MPa (20 bar), 1050 K: rho = P * M / (R_gas * T)
    # R_gas = 8.314 J/(mol*K), P = 2.0e6 Pa, T = 1050 K
    # rho = 2.0e6 * (M * 1e-3) / (8.314 * 1050) kg/m3 = 0.2291 * M [kg/m3] = 0.0002291 * M [g/cm3]
    if coolant_density is None:
        rho_calc = 0.0002291 * m_mol  # g/cm3
    else:
        rho_calc = coolant_density

    coolant.set_density("g/cm3", rho_calc)
    if he_frac > 0.0:
        coolant.add_nuclide("He4", he_frac, "ao")
    if xe_frac > 0.0:
        coolant.add_element("Xe", xe_frac, "ao")

    return materials, fuel, cladding, coolant, moderator, reflector, absorber, rho_calc

def run_coolant_composition_sweep(work_dir="scratch/study_coolant"):
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()
    os.chdir(work_dir)

    print("=" * 80)
    print("STUDY 1: COOLANT GAS COMPOSITION SWEEP (M = 4 to 40 g/mol)")
    print("=" * 80)

    # Core configuration at optimal H/R = 1.85 (R = 17.5 cm, H = 32.5 cm, or R = 18.0 cm, H = 33.3 cm)
    # Let's test on Opt-5 base first (R=16.0, H=43.5) and then on H/R=1.85 (R=17.7, H=35.5)
    compositions = [
        {"name": "He-Xe 72/28 (M=40)", "he": 0.72, "xe": 0.28},
        {"name": "He-Xe 80/20 (M=29.5)", "he": 0.80, "xe": 0.20},
        {"name": "He-Xe 85/15 (M=23.1)", "he": 0.85, "xe": 0.15},
        {"name": "He-Xe 90/10 (M=16.7)", "he": 0.90, "xe": 0.10},
        {"name": "Pure He (M=4.0)", "he": 1.00, "xe": 0.00},
    ]

    p_base = ReactorParameters(
        core_radius=17.70, active_core_height=35.40,  # H/R = 2.0 (near optimal buckling)
        moderator_fuel_pitch=1.58,                    # Slightly tighter pitch (V_mod/V_fuel = 0.98)
        radial_reflector_thickness=9.60,
        axial_reflector_thickness=7.60,
        control_drum_radius=4.30,
        absorber_thickness=1.45,
        absorber_sector_angular_width=140.0,
        b10_enrichment=0.90,
        particles=12000,
        batches=40,
        inactive=10
    )

    results = []
    # First: evaluate vacuum/voided state once
    print("\n--- Evaluating Voided Core State (Vacuum) ---")
    openmc.reset_auto_ids()
    mats_v, f_v, cl_v, co_v, mo_v, be_v, b_v, _ = build_custom_materials(p_base, he_frac=0.72, xe_frac=0.28, coolant_density=1.0e-7)
    geom_v, _, _ = build_geometry(p_base, f_v, cl_v, co_v, mo_v, be_v, b_v)
    sett_v = build_settings(p_base)
    mod_v = openmc.Model(materials=mats_v, geometry=geom_v, settings=sett_v)
    mod_v.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{p_base.batches}.h5") as sp:
        k_void = float(sp.keff.nominal_value)
        s_void = float(sp.keff.std_dev)
    print(f"Voided keff = {k_void:.5f} +/- {s_void:.5f}")

    for comp in compositions:
        openmc.reset_auto_ids()
        mats, f, cl, co, m, be, b, rho = build_custom_materials(p_base, he_frac=comp["he"], xe_frac=comp["xe"])
        geom, _, _ = build_geometry(p_base, f, cl, co, m, be, b)
        sett = build_settings(p_base)
        model = openmc.Model(materials=mats, geometry=geom, settings=sett)
        model.export_to_xml()
        openmc.run(output=False)
        with openmc.StatePoint(f"statepoint.{p_base.batches}.h5") as sp:
            k_op = float(sp.keff.nominal_value)
            s_op = float(sp.keff.std_dev)

        drho_void = ((k_void - 1.0)/k_void - (k_op - 1.0)/k_op) * 1e5
        u_void = np.sqrt((s_void/k_void**2)**2 + (s_op/k_op**2)**2) * 1e5
        conf_3sigma = drho_void + 3.0 * u_void

        rec = {
            "name": comp["name"], "he": comp["he"], "xe": comp["xe"],
            "density_g_cm3": rho, "k_op": k_op, "s_op": s_op,
            "k_void": k_void, "s_void": s_void,
            "drho_void_pcm": drho_void, "u_void_pcm": u_void,
            "conf_3sigma_pcm": conf_3sigma
        }
        results.append(rec)
        print(f"{comp['name']:24s} | rho = {rho:.5f} g/cm3 | k_op = {k_op:.5f} +/- {s_op:.5f} | Void drho = {drho_void:+.1f} +/- {u_void:.1f} pcm (3-sigma bound: {conf_3sigma:+.1f} pcm)")

    os.chdir(orig_dir)
    csv_path = "results/opt6/coolant_composition_void_study.csv"
    os.makedirs(os.path.dirname(csv_path), exist_ok=True)
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Composition", "He_ao", "Xe_ao", "Coolant Density [g/cm3]", "k_op", "s_op", "k_void", "s_void", "Void drho [pcm]", "u_void [pcm]", "3-sigma upper bound [pcm]"])
        for r in results:
            writer.writerow([
                r["name"], f"{r['he']:.2f}", f"{r['xe']:.2f}", f"{r['density_g_cm3']:.5f}",
                f"{r['k_op']:.5f}", f"{r['s_op']:.5f}", f"{r['k_void']:.5f}", f"{r['s_void']:.5f}",
                f"{r['drho_void_pcm']:+.1f}", f"{r['u_void_pcm']:.1f}", f"{r['conf_3sigma_pcm']:+.1f}"
            ])
    print(f"Saved: {csv_path}")
    return results

if __name__ == "__main__":
    run_coolant_composition_sweep()
