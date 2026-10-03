#!/usr/bin/env python3
"""
Targeted tuning for Opt-9 to achieve:
- k_out (BOL) ~ 1.018 to 1.025
- k_in (BOL) <= 0.955 (SDM > 4500 pcm)
- Delta-rho_void < -100 pcm
- Mass < 610 kg
"""

import os
import sys
import numpy as np
import openmc

os.environ["PATH"] = "/home/aleri/miniforge3/envs/openmc-env/bin:" + os.environ.get("PATH", "")
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
os.environ["OMP_NUM_THREADS"] = "8"
openmc.config['cross_sections'] = "/home/aleri/wmp_data/cross_sections_wmp.xml"

from run_opt8_master import Opt8Parameters, build_opt8_materials, build_opt8_geometry, build_opt8_settings

def test_geom(name, r_core, h_core, r_refl_rad=10.2, h_refl_ax=9.5, ch_outer=0.675, duct_r=0.72):
    p = Opt8Parameters()
    p.core_radius = r_core
    p.active_core_height = h_core
    p.radial_reflector_thickness = r_refl_rad
    p.axial_reflector_thickness = h_refl_ax
    p.radial_reflector_outer_radius = p.core_radius + p.radial_reflector_thickness
    p.total_reactor_height = p.active_core_height + 2 * p.axial_reflector_thickness
    p.channel_outer_radius = ch_outer
    p.drum_radial_position = p.core_radius + p.drum_clearance + p.control_drum_radius

    # Modify build_opt8_geometry inline to support configurable duct_r
    def build_geom(drum_angles, is_voided=False):
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
            cdk.rotation = (0.0, 0.0, k * 360.0 / p.number_of_control_drums + drum_angles[k])

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
    sett.particles = 50000
    sett.batches = 50
    sett.inactive = 12
    sett.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}

    # 1. k_out
    geom1, mats1, n_pins = build_geom([0.0]*8, is_voided=False)
    m1 = openmc.Model(materials=mats1, geometry=geom1, settings=sett)
    m1.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
        k_out = float(sp.keff.nominal_value)
        s_out = float(sp.keff.std_dev)

    # 2. k_void
    geom2, mats2, _ = build_geom([0.0]*8, is_voided=True)
    m2 = openmc.Model(materials=mats2, geometry=geom2, settings=sett)
    m2.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
        k_vd = float(sp.keff.nominal_value)
        s_vd = float(sp.keff.std_dev)

    # 3. k_in
    geom3, mats3, _ = build_geom([180.0]*8, is_voided=False)
    m3 = openmc.Model(materials=mats3, geometry=geom3, settings=sett)
    m3.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
        k_in = float(sp.keff.nominal_value)
        s_in = float(sp.keff.std_dev)

    drho_void = (1.0/k_out - 1.0/k_vd) * 1e5
    u_void = np.sqrt((s_out/k_out**2)**2 + (s_vd/k_vd**2)**2) * 1e5
    bound_3s_void = drho_void + 3.0 * u_void
    sdm = (1.0 - k_in) * 1e5
    worth = (1.0/k_in - 1.0/k_out) * 1e5

    print(f"\n{'='*70}")
    print(f"DESIGN: {name}")
    print(f"  R_core={r_core:.2f} cm, H_core={h_core:.2f} cm (H/R={h_core/r_core:.2f}), Pins={n_pins}")
    print(f"  Duct_r={duct_r:.3f} cm, Vol={vol_tot:.1f} L, Mass_est={mass_est:.1f} kg")
    print(f"  k_out (0°)   = {k_out:.5f} +/- {s_out:.5f}")
    print(f"  k_in  (180°) = {k_in:.5f} +/- {s_in:.5f}")
    print(f"  SDM          = {sdm:.0f} pcm | Drum Worth = {worth:.0f} pcm")
    print(f"  Void drho    = {drho_void:+.1f} +/- {u_void:.1f} pcm | 3-sigma bound = {bound_3s_void:+.1f} pcm")
    print(f"{'='*70}\n")
    return {
        "name": name, "k_out": k_out, "s_out": s_out, "k_in": k_in, "s_in": s_in,
        "drho_void": drho_void, "u_void": u_void, "bound_3s": bound_3s_void,
        "vol": vol_tot, "mass": mass_est, "pins": n_pins
    }

if __name__ == "__main__":
    # Test cases with slightly taller/wider cores to get k_out ~ 1.020:
    # Option 1: R=18.4 cm (385 pins), H=39.0 cm, duct_r=0.72 cm
    test_geom("Opt9-T1 (R18.4_H39.0_Duct0.72)", 18.40, 39.00, duct_r=0.72)
    # Option 2: R=18.4 cm (385 pins), H=40.0 cm, duct_r=0.75 cm
    test_geom("Opt9-T2 (R18.4_H40.0_Duct0.75)", 18.40, 40.00, duct_r=0.75)
