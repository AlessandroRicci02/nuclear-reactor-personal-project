"""
Test 10-Year Depletion Calculation for NEP Reactor
==================================================
Simulates continuous 10-year (3650 days) full-power depletion at 2.0 MWth.
"""
import os
import numpy as np
import openmc
import openmc.deplete

os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
from model import ReactorParameters, build_materials, build_geometry, build_settings, calculate_inventory
from test_drums import build_geometry_custom

def run_test_10yr_depletion(r_core=17.2, h_core=44.0, pitch=1.62, n_drums=8, r_drum=4.5, arc_deg=130.0, t_abs=1.2):
    p = ReactorParameters(
        core_radius=r_core,
        active_core_height=h_core,
        moderator_fuel_pitch=pitch,
        radial_reflector_thickness=10.5,
        control_drum_radius=r_drum,
        absorber_thickness=t_abs,
        absorber_sector_angular_width=arc_deg,
        particles=1500,
        batches=25,
        inactive=8,
        drum_angle_deg=0.0
    )
    p.number_of_control_drums = n_drums
    p.drum_radial_position = p.core_radius + 0.3 + p.control_drum_radius
    p.radial_reflector_outer_radius = p.drum_radial_position + p.control_drum_radius + 1.0
    p.radial_reflector_thickness = p.radial_reflector_outer_radius - p.core_radius

    openmc.reset_auto_ids()
    mats, fuel, clad, cool, mod, beo, b4c = build_materials(p)
    geom, pins, drums = build_geometry_custom(p, fuel, clad, cool, mod, beo, b4c)

    inv = calculate_inventory(p, len(pins))
    v_fuel_cm3 = inv["items"][0]["volume_cm3"]
    m_fuel_kg = inv["items"][0]["mass_kg"]
    fuel.volume = v_fuel_cm3
    fuel.depletable = True

    m_U_kg = m_fuel_kg * (237.457 / 251.464)
    m_u235_kg = m_U_kg * 0.1975

    power_W = 2.0e6
    settings = build_settings(p)
    model = openmc.Model(materials=mats, geometry=geom, settings=settings)

    chain_path = "chain_nep.xml"
    # Timesteps for 10 years: 30d, 150d, 185d (yr 1), 365d (yr 2), 730d (yr 4), 730d (yr 6), 730d (yr 8), 730d (yr 10) -> total 3650 days
    timesteps_d = [30.0, 150.0, 185.0, 365.0, 730.0, 730.0, 730.0, 730.0]

    print(f"Core: R={r_core}cm, H={h_core}cm, Pins={len(pins)}, U-235={m_u235_kg:.2f}kg")
    print(f"Starting 10-year depletion with {len(timesteps_d)} steps...")
    operator = openmc.deplete.CoupledOperator(model, chain_path)
    integrator = openmc.deplete.PredictorIntegrator(operator, timesteps_d, power=power_W, timestep_units='d')
    integrator.integrate()

    res = openmc.deplete.Results("depletion_results.h5")
    raw_time, k_vals = res.get_keff()
    times_days = np.array([t / 86400.0 for t in raw_time])
    keff_series = np.array([float(val[0]) for val in k_vals])
    keff_stds = np.array([float(val[1]) for val in k_vals])

    _, u235_g = res.get_mass(fuel, "U235")
    _, pu239_g = res.get_mass(fuel, "Pu239")

    print("\n10-Year Depletion History:")
    for t, k, s, u, pu in zip(times_days, keff_series, keff_stds, u235_g, pu239_g):
        print(f"  t = {t:6.1f} d ({t/365.25:5.2f} yr) | keff = {k:.5f} +/- {s:.5f} | U-235 = {u/1000:6.3f} kg | Pu-239 = {pu:6.1f} g")

    # Clean intermediate files
    for f in os.listdir("."):
        if f.startswith("openmc_simulation_n") and f.endswith(".h5"):
            try:
                os.remove(f)
            except:
                pass

if __name__ == "__main__":
    run_test_10yr_depletion()
