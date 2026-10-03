#!/usr/bin/env python3
"""
Test Candidate Opt-9 Geometry Variations to Find Optimal Balance:
- High negative void reactivity (Delta-rho_void < -100 pcm with 3-sigma bound < 0)
- Lifetime > 10 yr (k_out ~ 1.020 to 1.026)
- Shutdown and stuck-drum subcriticality
- Compact mass (target < 610 kg)
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

def run_candidate(name, r_core, h_core, r_refl_rad, h_refl_ax, duct_r):
    p = Opt8Parameters()
    p.core_radius = r_core
    p.active_core_height = h_core
    p.radial_reflector_thickness = r_refl_rad
    p.axial_reflector_thickness = h_refl_ax
    p.radial_reflector_outer_radius = p.core_radius + p.radial_reflector_thickness
    p.total_reactor_height = p.active_core_height + 2 * p.axial_reflector_thickness
    p.drum_radial_position = p.core_radius + p.drum_clearance + p.control_drum_radius

    # Materials
    openmc.reset_auto_ids()
    mats, fuel, cladding, coolant, moderator, reflector, absorber = build_opt8_materials(p)

    # Geometry builder with specified duct_r
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
        c_core = openmc.Cell(fill=lattice, region=-s_core & +z_refl_bot & -z_refl_top)

        # Control drums
        drum_cells = []
        drum_region_complement = []
        for d_idx in range(p.number_of_control_drums):
            phi = 2.0 * np.pi * d_idx / p.number_of_control_drums
            cx = p.drum_radial_position * np.cos(phi)
            cy = p.drum_radial_position * np.sin(phi)
            s_drum_cyl = openmc.ZCylinder(x0=cx, y0=cy, r=p.control_drum_radius)
            drum_region_complement.append(+s_drum_cyl)

            theta_deg = drum_angle
            theta_rad = np.radians(theta_deg)
            bisector_angle = phi + np.pi + theta_rad
            half_angle = np.radians(p.absorber_sector_angular_width / 2.0)
            ang1 = bisector_angle - half_angle
            ang2 = bisector_angle + half_angle

            p1 = openmc.Plane(a=np.sin(ang1), b=-np.cos(ang1), c=0.0, d=cx*np.sin(ang1) - cy*np.cos(ang1))
            p2 = openmc.Plane(a=np.sin(ang2), b=-np.cos(ang2), c=0.0, d=cx*np.sin(ang2) - cy*np.cos(ang2))
            s_abs_inner = openmc.ZCylinder(x0=cx, y0=cy, r=p.control_drum_radius - p.absorber_thickness)

            abs_sector = -s_drum_cyl & +s_abs_inner & -p1 & +p2
            c_abs = openmc.Cell(fill=absorber, region=abs_sector & +z_refl_bot & -z_refl_top)
            c_be = openmc.Cell(fill=reflector, region=(-s_drum_cyl & ~abs_sector) & +z_refl_bot & -z_refl_top)
            drum_cells.extend([c_abs, c_be])

        s_outer = openmc.ZCylinder(r=p.radial_reflector_outer_radius, boundary_type="vacuum")
        refl_region = +s_core & -s_outer & +z_refl_bot & -z_refl_top
        for comp in drum_region_complement:
            refl_region = refl_region & comp
        c_refl_matrix = openmc.Cell(fill=reflector, region=refl_region)

        geom = openmc.Geometry([c_core, c_refl_matrix] + drum_cells)
        return geom, mats, n_pins

    # Mass & volume calculation
    vol_tot = np.pi * p.radial_reflector_outer_radius**2 * p.total_reactor_height * 1e-3
    # More accurate mass:
    # Core vol:
    v_core = np.pi * p.core_radius**2 * p.active_core_height
    # Radial refl:
    v_refl_rad = np.pi * (p.radial_reflector_outer_radius**2 - p.core_radius**2) * p.total_reactor_height
    # Axial refl:
    v_refl_ax = np.pi * p.core_radius**2 * 2 * p.axial_reflector_thickness
    mass_est = (v_core * 5.5 + v_refl_rad * 2.85 + v_refl_ax * 2.85) * 1e-3 # kg approx

    sett = openmc.Settings()
    sett.particles = 40000
    sett.batches = 45
    sett.inactive = 10
    sett.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}

    # Run 1: Operating Drums OUT (0 deg)
    geom_out, mats_out, n_pins = build_geom(0.0, is_voided=False)
    m1 = openmc.Model(materials=mats_out, geometry=geom_out, settings=sett)
    m1.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
        k_out = float(sp.keff.nominal_value)
        s_out = float(sp.keff.std_dev)

    # Run 2: Voided Drums OUT
    geom_vd, mats_vd, _ = build_geom(0.0, is_voided=True)
    m2 = openmc.Model(materials=mats_vd, geometry=geom_vd, settings=sett)
    m2.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
        k_vd = float(sp.keff.nominal_value)
        s_vd = float(sp.keff.std_dev)

    # Run 3: Operating Drums IN (180 deg)
    geom_in, mats_in, _ = build_geom(180.0, is_voided=False)
    m3 = openmc.Model(materials=mats_in, geometry=geom_in, settings=sett)
    m3.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett.batches}.h5") as sp:
        k_in = float(sp.keff.nominal_value)
        s_in = float(sp.keff.std_dev)

    drho_void = (1.0/k_out - 1.0/k_vd) * 1e5
    u_void = np.sqrt((s_out/k_out**2)**2 + (s_vd/k_vd**2)**2) * 1e5
    bound_3s_void = drho_void + 3.0 * u_void
    drum_worth = (1.0/k_in - 1.0/k_out) * 1e5

    print(f"\n{'='*70}")
    print(f"CANDIDATE: {name}")
    print(f"  R_core={r_core:.2f} cm, H_core={h_core:.2f} cm (H/R={h_core/r_core:.2f}), Pins={n_pins}")
    print(f"  Refl_rad={r_refl_rad:.2f} cm, Refl_ax={h_refl_ax:.2f} cm, Duct_r={duct_r:.3f} cm")
    print(f"  Vol={vol_tot:.1f} L, Mass_est={mass_est:.1f} kg")
    print(f"  k_out = {k_out:.5f} +/- {s_out:.5f}")
    print(f"  k_in  = {k_in:.5f} +/- {s_in:.5f}")
    print(f"  Drum Worth = {drum_worth:.0f} pcm")
    print(f"  Void drho  = {drho_void:+.1f} +/- {u_void:.1f} pcm | 3-sigma bound = {bound_3s_void:+.1f} pcm")
    print(f"{'='*70}\n")
    return {
        "name": name, "r_core": r_core, "h_core": h_core, "pins": n_pins,
        "k_out": k_out, "s_out": s_out, "k_in": k_in, "s_in": s_in,
        "drho_void": drho_void, "u_void": u_void, "bound_3s": bound_3s_void,
        "vol": vol_tot, "mass": mass_est
    }

if __name__ == "__main__":
    candidates = [
        # Candidate A: R=18.6 cm (397 pins), H=36.0 cm, Duct=0.72 cm, Ax_Refl=9.5 cm, Rad_Refl=10.0 cm
        {"name": "Cand-A (R18.6_H36.0_Duct0.72)", "r_core": 18.60, "h_core": 36.00, "r_refl_rad": 10.00, "h_refl_ax": 9.50, "duct_r": 0.72},
        # Candidate B: R=18.6 cm (397 pins), H=37.5 cm, Duct=0.75 cm, Ax_Refl=9.5 cm, Rad_Refl=10.0 cm
        {"name": "Cand-B (R18.6_H37.5_Duct0.75)", "r_core": 18.60, "h_core": 37.50, "r_refl_rad": 10.00, "h_refl_ax": 9.50, "duct_r": 0.75},
        # Candidate C: R=18.2 cm (373 pins), H=38.0 cm, Duct=0.72 cm, Ax_Refl=9.5 cm, Rad_Refl=10.2 cm
        {"name": "Cand-C (R18.2_H38.0_Duct0.72)", "r_core": 18.20, "h_core": 38.00, "r_refl_rad": 10.20, "h_refl_ax": 9.50, "duct_r": 0.72},
    ]

    for c in candidates:
        run_candidate(c["name"], c["r_core"], c["h_core"], c["r_refl_rad"], c["h_refl_ax"], c["duct_r"])
