"""
Test Drum Parameters & Geometry
"""
import os
import numpy as np
import openmc

os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"

from model import ReactorParameters, build_materials, build_geometry, build_settings

def test_config(n_drums, r_drum, arc_deg, t_abs, clearance, r_core=17.2, h_core=44.0, pitch=1.62):
    p = ReactorParameters(
        core_radius=r_core,
        active_core_height=h_core,
        moderator_fuel_pitch=pitch,
        radial_reflector_thickness=10.5,
        control_drum_radius=r_drum,
        absorber_thickness=t_abs,
        absorber_sector_angular_width=arc_deg,
        particles=3000,
        batches=25,
        inactive=8
    )
    p.number_of_control_drums = n_drums
    p.drum_radial_position = p.core_radius + clearance + p.control_drum_radius
    p.radial_reflector_outer_radius = max(p.radial_reflector_outer_radius, p.drum_radial_position + p.control_drum_radius + 1.0)
    p.radial_reflector_thickness = p.radial_reflector_outer_radius - p.core_radius

    # 1. Drums out (0 deg)
    openmc.reset_auto_ids()
    p.drum_angle_deg = 0.0
    mats, fuel, clad, cool, mod, beo, b4c = build_materials(p)

    # Geometry with parameterized arc
    geom, pins, drums = build_geometry_custom(p, fuel, clad, cool, mod, beo, b4c)
    settings = build_settings(p)
    model = openmc.Model(materials=mats, geometry=geom, settings=settings)
    model.export_to_xml()
    openmc.run(output=False)

    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_out = float(sp.keff.nominal_value)
        k_out_std = float(sp.keff.std_dev)

    # 2. Drums in (180 deg)
    openmc.reset_auto_ids()
    p.drum_angle_deg = 180.0
    mats, fuel, clad, cool, mod, beo, b4c = build_materials(p)
    geom, pins, drums = build_geometry_custom(p, fuel, clad, cool, mod, beo, b4c)
    settings = build_settings(p)
    model = openmc.Model(materials=mats, geometry=geom, settings=settings)
    model.export_to_xml()
    openmc.run(output=False)

    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_in = float(sp.keff.nominal_value)
        k_in_std = float(sp.keff.std_dev)

    # 3. 1 Stuck Drum
    openmc.reset_auto_ids()
    stuck_angles = [0.0] + [180.0] * (n_drums - 1)
    mats, fuel, clad, cool, mod, beo, b4c = build_materials(p)
    geom, pins, drums = build_geometry_custom(p, fuel, clad, cool, mod, beo, b4c, drum_angles=stuck_angles)
    settings = build_settings(p)
    model = openmc.Model(materials=mats, geometry=geom, settings=settings)
    model.export_to_xml()
    openmc.run(output=False)

    with openmc.StatePoint(f"statepoint.{p.batches}.h5") as sp:
        k_stuck = float(sp.keff.nominal_value)
        k_stuck_std = float(sp.keff.std_dev)

    worth_pcm = (1.0 / k_in - 1.0 / k_out) * 1.0e5
    sdm_norm = (1.0 - k_in) / k_in * 1.0e5
    sdm_stk = (1.0 - k_stuck) / k_stuck * 1.0e5

    print(f"Results for N={n_drums}, R_drum={r_drum:.1f}, Arc={arc_deg:.0f}°, T_abs={t_abs:.2f}cm, Clear={clearance:.1f}cm:")
    print(f"  k_out = {k_out:.5f} +/- {k_out_std:.5f}")
    print(f"  k_in  = {k_in:.5f} +/- {k_in_std:.5f} (SDM = {sdm_norm:.0f} pcm)")
    print(f"  k_stk = {k_stuck:.5f} +/- {k_stuck_std:.5f} (SDM_stuck = {sdm_stk:.0f} pcm)")
    print(f"  Worth = {worth_pcm:.0f} pcm ({worth_pcm/1000:.2f}% Delta-k/k)\n")
    return k_out, k_in, k_stuck, worth_pcm

def build_geometry_custom(p, fuel, cladding, coolant, moderator, reflector, absorber, drum_angles=None):
    if drum_angles is None:
        drum_angles = [p.drum_angle_deg] * p.number_of_control_drums

    s_fuel = openmc.ZCylinder(r=p.fuel_radius, name="Pin Fuel Outer")
    s_clad = openmc.ZCylinder(r=p.clad_outer_radius, name="Pin Clad Outer")
    s_channel = openmc.ZCylinder(r=p.channel_outer_radius, name="Coolant Channel Outer")

    c_pin_fuel = openmc.Cell(name="Pin Fuel", fill=fuel, region=-s_fuel)
    c_pin_clad = openmc.Cell(name="Pin Clad", fill=cladding, region=+s_fuel & -s_clad)
    c_pin_cool = openmc.Cell(name="Pin Coolant", fill=coolant, region=+s_clad & -s_channel)
    c_pin_mod = openmc.Cell(name="Pin Mod", fill=moderator, region=+s_channel)
    u_pin = openmc.Universe(cells=[c_pin_fuel, c_pin_clad, c_pin_cool, c_pin_mod], name="Cooled Pin Universe")

    c_mod_only = openmc.Cell(name="Mod Only", fill=moderator)
    u_mod_only = openmc.Universe(cells=[c_mod_only], name="Moderator Universe")

    grid_dim = int(np.ceil(2 * p.core_radius / p.moderator_fuel_pitch))
    if grid_dim % 2 == 0:
        grid_dim += 1

    intact_pin_centers = []
    lattice_universes = []
    for i in range(grid_dim):
        row = []
        for j in range(grid_dim):
            x = -grid_dim * p.moderator_fuel_pitch / 2 + (j + 0.5) * p.moderator_fuel_pitch
            y = -grid_dim * p.moderator_fuel_pitch / 2 + (i + 0.5) * p.moderator_fuel_pitch
            if x**2 + y**2 <= (p.core_radius - p.channel_outer_radius)**2:
                row.append(u_pin)
                intact_pin_centers.append((x, y))
            else:
                row.append(u_mod_only)
        lattice_universes.append(row)

    lattice = openmc.RectLattice(name="Core Cooled Fuel Lattice")
    lattice.lower_left = (-grid_dim * p.moderator_fuel_pitch / 2, -grid_dim * p.moderator_fuel_pitch / 2)
    lattice.pitch = (p.moderator_fuel_pitch, p.moderator_fuel_pitch)
    lattice.universes = lattice_universes
    lattice.outer = u_mod_only

    # Control Drum with variable sector
    r_abs_in = p.control_drum_radius - p.absorber_thickness
    drum_outer_cyl = openmc.ZCylinder(r=p.control_drum_radius)
    drum_inner_cyl = openmc.ZCylinder(r=r_abs_in)

    half_phi_rad = np.radians(p.absorber_sector_angular_width / 2.0)
    sin_phi = np.sin(half_phi_rad)
    cos_phi = np.cos(half_phi_rad)
    p_sec1 = openmc.Plane(a=sin_phi, b=-cos_phi, c=0.0, d=0.0)
    p_sec2 = openmc.Plane(a=sin_phi, b=+cos_phi, c=0.0, d=0.0)
    p_x_pos = openmc.XPlane(x0=0.0)

    abs_region = -drum_outer_cyl & +drum_inner_cyl & +p_sec1 & +p_sec2 & +p_x_pos
    beo_region = -drum_outer_cyl & ~abs_region

    c_drum_abs = openmc.Cell(name="Drum B4C Sector", fill=absorber, region=abs_region)
    c_drum_beo = openmc.Cell(name="Drum BeO Body", fill=reflector, region=beo_region)
    u_drum = openmc.Universe(cells=[c_drum_abs, c_drum_beo], name="Control Drum Universe")

    s_core = openmc.ZCylinder(r=p.core_radius)
    s_outer = openmc.ZCylinder(r=p.radial_reflector_outer_radius, boundary_type="vacuum")

    z_core_bot = openmc.ZPlane(z0=-p.active_core_height / 2)
    z_core_top = openmc.ZPlane(z0=p.active_core_height / 2)
    z_refl_bot = openmc.ZPlane(z0=-p.total_reactor_height / 2, boundary_type="vacuum")
    z_refl_top = openmc.ZPlane(z0=p.total_reactor_height / 2, boundary_type="vacuum")

    drum_cells = []
    drum_regions = []
    drum_centers = []

    for k in range(p.number_of_control_drums):
        angle_rad = k * 2.0 * np.pi / p.number_of_control_drums
        drum_nominal_deg = k * 360.0 / p.number_of_control_drums
        xk = p.drum_radial_position * np.cos(angle_rad)
        yk = p.drum_radial_position * np.sin(angle_rad)
        drum_centers.append((xk, yk, drum_nominal_deg, drum_angles[k]))

        total_rot_deg = drum_nominal_deg + drum_angles[k]
        dk = openmc.ZCylinder(x0=xk, y0=yk, r=p.control_drum_radius)
        cdk = openmc.Cell(name=f"Control Drum {k}", fill=u_drum,
                          region=-dk & +z_core_bot & -z_core_top)
        cdk.translation = (xk, yk, 0.0)
        cdk.rotation = (0.0, 0.0, total_rot_deg)

        drum_cells.append(cdk)
        drum_regions.append(-dk)

    drum_union = drum_regions[0]
    for reg in drum_regions[1:]:
        drum_union = drum_union | reg

    c_core = openmc.Cell(name="Active Core Cell", fill=lattice,
                         region=-s_core & +z_core_bot & -z_core_top)
    c_rad_refl = openmc.Cell(name="Radial Reflector Matrix", fill=reflector,
                             region=+s_core & -s_outer & +z_core_bot & -z_core_top & ~drum_union)
    c_bot_refl = openmc.Cell(name="Bottom Axial Reflector", fill=reflector,
                             region=-s_outer & +z_refl_bot & -z_core_bot)
    c_top_refl = openmc.Cell(name="Top Axial Reflector", fill=reflector,
                             region=-s_outer & +z_core_top & -z_refl_top)

    root_universe = openmc.Universe(
        name="Root Reactor Universe",
        cells=[c_core, c_rad_refl, c_bot_refl, c_top_refl] + drum_cells
    )
    return openmc.Geometry(root_universe), intact_pin_centers, drum_centers

if __name__ == "__main__":
    # Test 1: 6 drums, R_drum=4.8, Arc=130, T_abs=1.2, Clear=0.3
    test_config(n_drums=6, r_drum=4.8, arc_deg=130.0, t_abs=1.2, clearance=0.3)

    # Test 2: 8 drums, R_drum=4.5, Arc=130, T_abs=1.2, Clear=0.3
    test_config(n_drums=8, r_drum=4.5, arc_deg=130.0, t_abs=1.2, clearance=0.3)
