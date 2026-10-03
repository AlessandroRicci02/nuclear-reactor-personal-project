#!/usr/bin/env python3
"""
Test wider coolant channels and larger BeO streaming ducts:
- Clad outer: 0.63 cm
- Channel outer: 0.70 cm (coolant gap 0.70 mm vs 0.45 mm)
- Duct radius in BeO: 0.80 cm
- Test R_core = 18.2 to 18.6 cm, H_core = 38.0 to 40.0 cm
"""

import os
import sys
import numpy as np
import openmc

os.environ["PATH"] = "/home/aleri/miniforge3/envs/openmc-env/bin:" + os.environ.get("PATH", "")
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
os.environ["OMP_NUM_THREADS"] = "8"
openmc.config['cross_sections'] = "/home/aleri/wmp_data/cross_sections_wmp.xml"

from run_opt8_master import Opt8Parameters, build_opt8_materials

def run_streaming_test(name, r_core, h_core, ch_outer, duct_r, p_MPa=2.5):
    p = Opt8Parameters()
    p.core_radius = r_core
    p.active_core_height = h_core
    p.channel_outer_radius = ch_outer
    p.radial_reflector_outer_radius = p.core_radius + p.radial_reflector_thickness
    p.total_reactor_height = p.active_core_height + 2 * p.axial_reflector_thickness
    p.drum_radial_position = p.core_radius + p.drum_clearance + p.control_drum_radius
    p.coolant_pressure_MPa = p_MPa
    p.coolant_density = (p.coolant_pressure_MPa * 1.0e6 * 4.0026e-3) / (8.31446 * p.coolant_temp_K) * 1.0e-3

    def build_geom(drum_angle, is_voided=False):
        openmc.reset_auto_ids()
        mats, fuel, cladding, coolant, moderator, reflector, absorber = build_opt8_materials(p)
        if is_voided:
            coolant.set_density("g/cm3", 1.0e-7)

        s_fuel = openmc.ZCylinder(r=p.fuel_radius)
        s_clad = openmc.ZCylinder(r=p.clad_outer_radius)
        s_channel = openmc.ZCylinder(r=p.channel_outer_radius)

        z_core_bot = openmc.ZPlane(z0=-p.active_core_height / 2.0)
        z_core_top = openmc.ZPlane(z0=+p.active_core_height / 2.0)
        z_refl_bot = openmc.ZPlane(z0=-p.total_reactor_height / 2.0, boundary_type="vacuum")
        z_refl_top = openmc.ZPlane(z0=+p.total_reactor_height / 2.0, boundary_type="vacuum")

        c_f = openmc.Cell(fill=fuel, region=-s_fuel & +z_core_bot & -z_core_top)
        c_cl = openmc.Cell(fill=cladding, region=+s_fuel & -s_clad & +z_core_bot & -z_core_top)
        c_co = openmc.Cell(fill=coolant, region=+s_clad & -s_channel & +z_core_bot & -z_core_top)
        c_mo = openmc.Cell(fill=moderator, region=+s_channel & +z_core_bot & -z_core_top)

        s_duct = openmc.ZCylinder(r=duct_r)
        c_co_bot = openmc.Cell(fill=coolant, region=-s_duct & +z_refl_bot & -z_core_bot)
        c_be_bot = openmc.Cell(fill=reflector, region=+s_duct & +z_refl_bot & -z_core_bot)
        c_co_top = openmc.Cell(fill=coolant, region=-s_duct & +z_core_top & -z_refl_top)
        c_be_top = openmc.Cell(fill=reflector, region=+s_duct & +z_core_top & -z_refl_top)

        u_pin = openmc.Universe(cells=[c_f, c_cl, c_co, c_mo, c_co_bot, c_be_bot, c_co_top, c_be_top])

        c_mo_core = openmc.Cell(fill=moderator, region=+z_core_bot & -z_core_top)
        c_be_bot_only = openmc.Cell(fill=reflector, region=+z_refl_bot & -z_core_bot)
        c_be_top_only = openmc.Cell(fill=reflector, region=+z_core_top & -z_refl_top)
        u_mod_only = openmc.Universe(cells=[c_mo_core, c_be_bot_only, c_be_top_only])

        grid_dim = int(np.ceil(2 * p.core_radius / p.moderator_fuel_pitch))
        if grid_dim % 2 == 0: grid_dim += 1
        lat_univs = []
        n_pins = 0
        for i in range(grid_dim):
            row = []
            for j in range(grid_dim):
                x = -grid_dim * p.moderator_fuel_pitch / 2.0 + (j + 0.5) * p.moderator_fuel_pitch
                y = -grid_dim * p.moderator_fuel_pitch / 2.0 + (i + 0.5) * p.moderator_fuel_pitch
                if x**2 + y**2 <= (p.core_radius - p.channel_outer_radius)**2:
                    row.append(u_pin)
                    n_pins += 1
                else:
                    row.append(u_mod_only)
            lat_univs.append(row)

        lattice = openmc.RectLattice()
        lattice.lower_left = (-grid_dim * p.moderator_fuel_pitch / 2.0, -grid_dim * p.moderator_fuel_pitch / 2.0)
        lattice.pitch = (p.moderator_fuel_pitch, p.moderator_fuel_pitch)
        lattice.universes = lat_univs
        lattice.outer = u_mod_only

        s_core = openmc.ZCylinder(r=p.core_radius)
        s_outer = openmc.ZCylinder(r=p.radial_reflector_outer_radius, boundary_type="vacuum")

        r_abs_in = p.control_drum_radius - p.absorber_thickness
        drum_outer_cyl = openmc.ZCylinder(r=p.control_drum_radius)
        drum_inner_cyl = openmc.ZCylinder(r=r_abs_in)

        half_phi = np.radians(p.absorber_sector_angular_width / 2.0)
        p1 = openmc.Plane(a=np.sin(half_phi), b=-np.cos(half_phi), c=0.0, d=0.0)
        p2 = openmc.Plane(a=np.sin(half_phi), b=+np.cos(half_phi), c=0.0, d=0.0)
        px = openmc.XPlane(x0=0.0)

        abs_region = -drum_outer_cyl & +drum_inner_cyl & +p1 & +p2 & +px
        beo_region = -drum_outer_cyl & ~abs_region

        c_abs = openmc.Cell(fill=absorber, region=abs_region)
        c_be = openmc.Cell(fill=reflector, region=beo_region)
        u_drum = openmc.Universe(cells=[c_abs, c_be])

        drum_cells = []
        drum_regions = []
        for k in range(p.number_of_control_drums):
            ang = k * 2.0 * np.pi / p.number_of_control_drums
            xk = p.drum_radial_position * np.cos(ang)
            yk = p.drum_radial_position * np.sin(ang)
            dk = openmc.ZCylinder(x0=xk, y0=yk, r=p.control_drum_radius)

            cdk = openmc.Cell(fill=u_drum, region=-dk & +z_core_bot & -z_core_top)
            cdk.translation = (xk, yk, 0.0)
            cdk.rotation = (0.0, 0.0, k * 360.0 / p.number_of_control_drums + drum_angle)

            drum_cells.append(cdk)
            drum_regions.append(-dk)

        drum_union = drum_regions[0]
        for reg in drum_regions[1:]:
            drum_union = drum_union | reg

        c_core_cyl = openmc.Cell(fill=lattice, region=-s_core & +z_refl_bot & -z_refl_top)
        c_rad_refl = openmc.Cell(fill=reflector, region=+s_core & -s_outer & +z_core_bot & -z_core_top & ~drum_union)
        c_rad_bot = openmc.Cell(fill=reflector, region=+s_core & -s_outer & +z_refl_bot & -z_core_bot)
        c_rad_top = openmc.Cell(fill=reflector, region=+s_core & -s_outer & +z_core_top & -z_refl_top)

        root_u = openmc.Universe(cells=[c_core_cyl, c_rad_refl, c_rad_bot, c_rad_top] + drum_cells)
        geom = openmc.Geometry(root_u)
        return geom, mats, n_pins

    vol_tot = np.pi * p.radial_reflector_outer_radius**2 * p.total_reactor_height * 1e-3
    mass_est = 628.49 * (vol_tot / 154.06)

    sett = openmc.Settings()
    sett.particles = 60000
    sett.batches = 55
    sett.inactive = 15
    sett.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}

    # 1. Operating Drums OUT
    geom_op, mats_op, n_pins = build_geom(0.0, is_voided=False)
    m_op = openmc.Model(materials=mats_op, geometry=geom_op, settings=sett)
    m_op.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
        k_op = float(sp.keff.nominal_value)
        s_op = float(sp.keff.std_dev)

    # 2. Voided Drums OUT
    geom_vd, mats_vd, _ = build_geom(0.0, is_voided=True)
    m_vd = openmc.Model(materials=mats_vd, geometry=geom_vd, settings=sett)
    m_vd.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
        k_vd = float(sp.keff.nominal_value)
        s_vd = float(sp.keff.std_dev)

    drho_void = (1.0/k_op - 1.0/k_vd) * 1e5
    u_void = np.sqrt((s_op/k_op**2)**2 + (s_vd/k_vd**2)**2) * 1e5
    bound_3s = drho_void + 3.0 * u_void

    print(f"\n======================================================================")
    print(f"CASE: {name}")
    print(f"  R_core={r_core:.2f} cm, H_core={h_core:.2f} cm, Pins={n_pins}")
    print(f"  Channel_r={ch_outer:.3f} cm, Duct_r={duct_r:.3f} cm, P={p_MPa:.1f} MPa")
    print(f"  Vol={vol_tot:.1f} L, Mass_est={mass_est:.1f} kg")
    print(f"  k_op = {k_op:.5f} +/- {s_op:.5f}")
    print(f"  k_vd = {k_vd:.5f} +/- {s_vd:.5f}")
    print(f"  Delta-rho_void = {drho_void:+.1f} +/- {u_void:.1f} pcm | 3-sigma bound: {bound_3s:+.1f} pcm")
    print(f"======================================================================\n")

if __name__ == "__main__":
    # Case 1: R=18.4 cm, H=41.0 cm, duct_r=0.80 cm, gap=0.70 mm (ch_outer=0.70 cm)
    run_streaming_test("WiderChannel_Duct0.80", 18.40, 41.00, ch_outer=0.700, duct_r=0.800)
    # Case 2: R=18.2 cm, H=41.8 cm, duct_r=0.85 cm, gap=0.55 mm (ch_outer=0.685 cm)
    run_streaming_test("Opt8Mod_Duct0.85", 18.20, 41.80, ch_outer=0.685, duct_r=0.850)
