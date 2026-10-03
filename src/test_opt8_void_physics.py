import openmc
import numpy as np
import sys
import csv
sys.path.append('src')
from model import ReactorParameters, build_materials, build_geometry, build_settings

cases = [
    {"name": "Pure He, 2.0 MPa, gap 0.45mm, P 1.612", "he": 1.0, "xe": 0.0, "P": 2.0, "gap": 0.045, "pitch": 1.612},
    {"name": "Pure He, 3.0 MPa, gap 0.45mm, P 1.612", "he": 1.0, "xe": 0.0, "P": 3.0, "gap": 0.045, "pitch": 1.612},
    {"name": "Pure He, 3.0 MPa, gap 0.50mm, P 1.600", "he": 1.0, "xe": 0.0, "P": 3.0, "gap": 0.050, "pitch": 1.600},
    {"name": "He-Xe 98/2, 2.5 MPa, gap 0.50mm, P 1.600", "he": 0.98, "xe": 0.02, "P": 2.5, "gap": 0.050, "pitch": 1.600},
]

results = []

for c in cases:
    p = ReactorParameters(
        core_radius=17.80, active_core_height=34.80,
        fuel_radius=0.58, cladding_thickness=0.05,
        coolant_channel_thickness=c["gap"], moderator_fuel_pitch=c["pitch"],
        radial_reflector_thickness=9.60, axial_reflector_thickness=7.50,
        number_of_control_drums=8, control_drum_radius=4.30,
        absorber_thickness=1.45, absorber_sector_angular_width=140.0,
        drum_clearance=0.20, b10_enrichment=0.90, drum_angle_deg=0.0,
        particles=60000, batches=60, inactive=15
    )

    M_cool = c["he"] * 4.0026 + c["xe"] * 131.293
    rho_cool = (c["P"] * 1.0e6 * M_cool * 1.0e-3) / (8.31446 * 1050.0) * 1.0e-3

    # Operating
    openmc.reset_auto_ids()
    m_op, f_op, cl_op, co_op, mo_op, be_op, b_op = build_materials(p)
    co_op.nuclides.clear()
    co_op.set_density('g/cm3', rho_cool)
    co_op.add_nuclide('He4', c["he"], 'ao')
    if c["xe"] > 0:
        co_op.add_element('Xe', c["xe"], 'ao')
    g_op, _, _ = build_geometry(p, f_op, cl_op, co_op, mo_op, be_op, b_op)
    mod_op = openmc.Model(materials=m_op, geometry=g_op, settings=build_settings(p))
    mod_op.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f'statepoint.{p.batches}.h5') as sp:
        k_op = float(sp.keff.nominal_value)
        s_op = float(sp.keff.std_dev)

    # Voided
    openmc.reset_auto_ids()
    m_vd, f_vd, cl_vd, co_vd, mo_vd, be_vd, b_vd = build_materials(p)
    co_vd.nuclides.clear()
    co_vd.set_density('g/cm3', 1.0e-7)
    co_vd.add_nuclide('He4', 1.0, 'ao')
    g_vd, _, _ = build_geometry(p, f_vd, cl_vd, co_vd, mo_vd, be_vd, b_vd)
    mod_vd = openmc.Model(materials=m_vd, geometry=g_vd, settings=build_settings(p))
    mod_vd.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f'statepoint.{p.batches}.h5') as sp:
        k_vd = float(sp.keff.nominal_value)
        s_vd = float(sp.keff.std_dev)

    drho = (1.0 / k_op - 1.0 / k_vd) * 1e5
    u_drho = np.sqrt((s_op / k_op**2)**2 + (s_vd / k_vd**2)**2) * 1e5
    bound_3sigma = drho + 3 * u_drho

    print(f"CASE: {c['name']}")
    print(f"  k_op = {k_op:.5f} +/- {s_op:.5f}")
    print(f"  k_vd = {k_vd:.5f} +/- {s_vd:.5f}")
    print(f"  Delta-rho_void = {drho:+.1f} +/- {u_drho:.1f} pcm")
    print(f"  3-sigma upper bound = {bound_3sigma:+.1f} pcm")

