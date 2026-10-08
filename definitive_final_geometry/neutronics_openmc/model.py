#!/usr/bin/env python3
"""
model.py - Parametric OpenMC CSG Model Generator for NEP Reactor Candidates
Supports dual-pass matrix cores with auxiliary channels and dual-annulus sleeve cores.
Includes temperature assignments, S(alpha,beta) scattering, and rotating control drums.
"""

from pathlib import Path
import math
import numpy as np
import openmc

XS_PATH = Path('/home/aleri/OpenMC_projects/first_test/new_possible_best_configurations_failed/campaign/data/endfb-viii.0-hdf5/cross_sections.xml')
OPENMC_BIN = '/home/aleri/miniforge3/envs/openmc-env/bin/openmc'

def get_hex_positions(n_rings, pitch_cm):
    """Generate (q, r, x, y, ring) for hexagonal lattice rings 0 to n_rings."""
    positions = []
    for q in range(-n_rings, n_rings + 1):
        for r in range(-n_rings, n_rings + 1):
            ring = max(abs(q), abs(r), abs(q + r))
            if ring <= n_rings:
                x = pitch_cm * (q + r / 2.0)
                y = pitch_cm * (math.sqrt(3.0) / 2.0 * r)
                positions.append((q, r, x, y, ring))
    # Sort by ring, then q, then r
    positions.sort(key=lambda item: (item[4], item[0], item[1]))
    return positions

def build_model(spec, theta_deg=0.0, T_fuel=1250.0, T_mod=900.0, T_cool=1000.0,
                T_refl=900.0, mod_density_factor=1.0, fuel_density_factor=1.0,
                f_in_mat=None, f_out_mat=None,
                particles=8000, batches=70, inactive=25, seed=42):
    """
    Build OpenMC Model for candidate specification.
    theta_deg: control drum rotation angle (0 = OUT, 180 = IN)
    mod_density_factor: multiplier on moderator density (for thermal expansion feedback)
    fuel_density_factor: multiplier on fuel density (for prompt thermal expansion feedback)
    f_in_mat, f_out_mat: optional burned fuel materials for depletion analysis
    """
    openmc.reset_auto_ids()
    openmc.config['cross_sections'] = str(XS_PATH)

    # -------------------------------------------------------------------------
    # 1. MATERIALS
    # -------------------------------------------------------------------------
    v_in = 19.0 * math.pi * (spec['r_fuel_cm']**2) * spec['height_cm']
    v_out = 198.0 * math.pi * (spec['r_fuel_cm']**2) * spec['height_cm']

    # Fuel UN Inner
    if f_in_mat is not None:
        mat_fuel_in = f_in_mat.clone()
        mat_fuel_in.id = 1
        mat_fuel_in.name = "UN_Inner"
        mat_fuel_in.temperature = T_fuel
    else:
        e_in = spec['enrichment_inner']
        a235_in = (e_in / 235.0439) / (e_in / 235.0439 + (1.0 - e_in) / 238.0508)
        a238_in = 1.0 - a235_in
        mat_fuel_in = openmc.Material(material_id=1, name="UN_Inner")
        mat_fuel_in.set_density("g/cm3", 13.50 * fuel_density_factor)
        mat_fuel_in.add_nuclide("U235", a235_in, "ao")
        mat_fuel_in.add_nuclide("U238", a238_in, "ao")
        mat_fuel_in.add_element("N", 1.0, "ao")
        mat_fuel_in.temperature = T_fuel
    mat_fuel_in.depletable = True
    mat_fuel_in.volume = v_in

    # Fuel UN Outer
    if f_out_mat is not None:
        mat_fuel_out = f_out_mat.clone()
        mat_fuel_out.id = 2
        mat_fuel_out.name = "UN_Outer"
        mat_fuel_out.temperature = T_fuel
    else:
        e_out = spec['enrichment_outer']
        a235_out = (e_out / 235.0439) / (e_out / 235.0439 + (1.0 - e_out) / 238.0508)
        a238_out = 1.0 - a235_out
        mat_fuel_out = openmc.Material(material_id=2, name="UN_Outer")
        mat_fuel_out.set_density("g/cm3", 13.50 * fuel_density_factor)
        mat_fuel_out.add_nuclide("U235", a235_out, "ao")
        mat_fuel_out.add_nuclide("U238", a238_out, "ao")
        mat_fuel_out.add_element("N", 1.0, "ao")
        mat_fuel_out.temperature = T_fuel
    mat_fuel_out.depletable = True
    mat_fuel_out.volume = v_out

    # Cladding Nb-1Zr
    mat_clad = openmc.Material(material_id=3, name="Nb1Zr")
    mat_clad.set_density("g/cm3", 8.57)
    mat_clad.add_element("Nb", 0.99, "wo")
    mat_clad.add_element("Zr", 0.01, "wo")
    mat_clad.temperature = (T_fuel + T_cool) / 2.0

    # Helium Coolant (Ideal gas at 2.5 MPa)
    P_Pa = spec['pressure_MPa'] * 1e6
    R_gas = 8314.46 / 4.0026 # J/(kg K)
    rho_he = (P_Pa / (R_gas * T_cool)) * 1e-3 # g/cm3
    mat_cool = openmc.Material(material_id=4, name="Helium")
    mat_cool.set_density("g/cm3", rho_he)
    mat_cool.add_element("He", 1.0, "ao")
    mat_cool.temperature = T_cool

    # Al2O3 Hydrogen Barrier / Ceramic Liner
    mat_al2o3 = openmc.Material(material_id=5, name="Al2O3")
    mat_al2o3.set_density("g/cm3", 3.95)
    mat_al2o3.add_element("Al", 2.0, "ao")
    mat_al2o3.add_element("O", 3.0, "ao")
    mat_al2o3.temperature = T_mod

    # YHx Moderator
    y_x = spec.get('y_stoich', 2.0)
    rho_mod = spec.get('rho_mod', 4.30) * mod_density_factor
    mat_mod = openmc.Material(material_id=6, name=f"YH_{y_x}")
    mat_mod.set_density("g/cm3", rho_mod)
    mat_mod.add_element("Y", 1.0, "ao")
    mat_mod.add_nuclide("H1", y_x, "ao")
    mat_mod.add_s_alpha_beta("c_H_in_YH2")
    mat_mod.add_s_alpha_beta("c_Y_in_YH2")
    mat_mod.temperature = T_mod

    # Reflector Material (BeO or Be)
    refl_type = spec.get('reflector', 'BeO')
    if refl_type == "BeO":
        mat_refl = openmc.Material(material_id=7, name="BeO")
        mat_refl.set_density("g/cm3", 2.90)
        mat_refl.add_element("Be", 1.0, "ao")
        mat_refl.add_element("O", 1.0, "ao")
        mat_refl.add_s_alpha_beta("c_Be")
    else:
        mat_refl = openmc.Material(material_id=7, name="Be")
        mat_refl.set_density("g/cm3", 1.85)
        mat_refl.add_element("Be", 1.0, "ao")
        mat_refl.add_s_alpha_beta("c_Be")
    mat_refl.temperature = T_refl

    # B4C Absorber (Enriched to 90 wt% B-10)
    a10 = (0.90 / 10.0129) / (0.90 / 10.0129 + 0.10 / 11.0093)
    a11 = 1.0 - a10
    mat_b4c = openmc.Material(material_id=8, name="B4C_Absorber")
    mat_b4c.set_density("g/cm3", 2.52)
    mat_b4c.add_nuclide("B10", 4.0 * a10, "ao")
    mat_b4c.add_nuclide("B11", 4.0 * a11, "ao")
    mat_b4c.add_element("C", 1.0, "ao")
    mat_b4c.temperature = T_refl

    materials = openmc.Materials([mat_fuel_in, mat_fuel_out, mat_clad, mat_cool,
                                  mat_al2o3, mat_mod, mat_refl, mat_b4c])

    # -------------------------------------------------------------------------
    # 2. GEOMETRY
    # -------------------------------------------------------------------------
    H_core = spec['height_cm']
    H_ax = spec['t_refl_ax_cm']
    z_bottom_core = -H_core / 2.0
    z_top_core = H_core / 2.0
    z_bottom_total = z_bottom_core - H_ax
    z_top_total = z_top_core + H_ax

    surf_z_bot_tot = openmc.ZPlane(z0=z_bottom_total, boundary_type='vacuum')
    surf_z_bot_core = openmc.ZPlane(z0=z_bottom_core)
    surf_z_top_core = openmc.ZPlane(z0=z_top_core)
    surf_z_top_tot = openmc.ZPlane(z0=z_top_total, boundary_type='vacuum')

    # Radial dimensions
    pitch = spec['pitch_cm']
    r_fuel = spec['r_fuel_cm']
    t_clad = spec['t_clad_cm']
    r_clad = r_fuel + t_clad
    t_hot_gap = spec['t_hot_gap_cm']
    r_hot_chan = r_clad + t_hot_gap
    t_al2o3 = spec.get('t_al2o3_cm', 0.010)
    r_liner = r_hot_chan + t_al2o3

    family = spec.get('family', 'matrix_dual_pass')

    # Calculate Core Radius avoiding pin clipping
    n_pins = spec['pins']
    n_rings = spec['rings']
    R_core = (n_rings + 0.57) * pitch
    drum_clearance = 0.18 # cm
    r_drum = spec['r_drum_cm']
    R_drum_center = R_core + drum_clearance + r_drum
    t_refl_rad = spec['t_refl_rad_cm']
    R_ref = R_core + t_refl_rad

    surf_cyl_core = openmc.ZCylinder(r=R_core)
    surf_cyl_ref = openmc.ZCylinder(r=R_ref, boundary_type='vacuum')

    # Surfaces for single pin concentric cylinders
    surf_fuel = openmc.ZCylinder(r=r_fuel)
    surf_clad = openmc.ZCylinder(r=r_clad)
    surf_hot = openmc.ZCylinder(r=r_hot_chan)
    surf_liner = openmc.ZCylinder(r=r_liner)

    # Unit cells: Inner UN vs Outer UN
    if family != "dual_annulus":
        aux_holes = spec.get('aux_holes', 6)
        r_aux = spec.get('r_aux_cm', 0.165)
        r_aux_pos = 0.40 * pitch # radial position of aux holes inside cell

        aux_cyls = []
        for a in range(aux_holes):
            ang = a * (2.0 * math.pi / aux_holes)
            ax_pos = r_aux_pos * math.cos(ang)
            ay_pos = r_aux_pos * math.sin(ang)
            aux_cyls.append(openmc.ZCylinder(x0=ax_pos, y0=ay_pos, r=r_aux))

        # Universe for Inner Pin Cell
        u_pin_in = openmc.Universe(name="Pin_Cell_Inner")
        c_f_in = openmc.Cell(fill=mat_fuel_in, region=-surf_fuel)
        c_cl_in = openmc.Cell(fill=mat_clad, region=+surf_fuel & -surf_clad)
        c_he_in = openmc.Cell(fill=mat_cool, region=+surf_clad & -surf_hot)
        c_al_in = openmc.Cell(fill=mat_al2o3, region=+surf_hot & -surf_liner)
        
        aux_cells_in = []
        mod_region_in = +surf_liner
        for ac in aux_cyls:
            aux_cells_in.append(openmc.Cell(fill=mat_cool, region=-ac))
            mod_region_in = mod_region_in & +ac
        c_mod_in = openmc.Cell(fill=mat_mod, region=mod_region_in)
        u_pin_in.add_cells([c_f_in, c_cl_in, c_he_in, c_al_in, c_mod_in] + aux_cells_in)

        # Universe for Outer Pin Cell
        u_pin_out = openmc.Universe(name="Pin_Cell_Outer")
        c_f_out = openmc.Cell(fill=mat_fuel_out, region=-surf_fuel)
        c_cl_out = openmc.Cell(fill=mat_clad, region=+surf_fuel & -surf_clad)
        c_he_out = openmc.Cell(fill=mat_cool, region=+surf_clad & -surf_hot)
        c_al_out = openmc.Cell(fill=mat_al2o3, region=+surf_hot & -surf_liner)
        
        aux_cells_out = []
        mod_region_out = +surf_liner
        for ac in aux_cyls:
            aux_cells_out.append(openmc.Cell(fill=mat_cool, region=-ac))
            mod_region_out = mod_region_out & +ac
        c_mod_out = openmc.Cell(fill=mat_mod, region=mod_region_out)
        u_pin_out.add_cells([c_f_out, c_cl_out, c_he_out, c_al_out, c_mod_out] + aux_cells_out)

    else:
        # Dual-annulus family
        t_part = spec.get('t_part_cm', 0.035)
        t_cold = spec.get('t_cold_gap_cm', 0.080)
        r_part_out = r_hot_chan + t_part
        r_cold_out = r_part_out + t_cold
        r_liner_da = r_cold_out + t_al2o3

        surf_part = openmc.ZCylinder(r=r_part_out)
        surf_cold = openmc.ZCylinder(r=r_cold_out)
        surf_liner_da = openmc.ZCylinder(r=r_liner_da)

        u_pin_in = openmc.Universe(name="Pin_Cell_Inner_DA")
        c_f = openmc.Cell(fill=mat_fuel_in, region=-surf_fuel)
        c_cl = openmc.Cell(fill=mat_clad, region=+surf_fuel & -surf_clad)
        c_he_hot = openmc.Cell(fill=mat_cool, region=+surf_clad & -surf_hot)
        c_part = openmc.Cell(fill=mat_clad, region=+surf_hot & -surf_part)
        c_he_cold = openmc.Cell(fill=mat_cool, region=+surf_part & -surf_cold)
        c_al = openmc.Cell(fill=mat_al2o3, region=+surf_cold & -surf_liner_da)
        c_mod = openmc.Cell(fill=mat_mod, region=+surf_liner_da)
        u_pin_in.add_cells([c_f, c_cl, c_he_hot, c_part, c_he_cold, c_al, c_mod])

        u_pin_out = openmc.Universe(name="Pin_Cell_Outer_DA")
        c_f2 = openmc.Cell(fill=mat_fuel_out, region=-surf_fuel)
        c_cl2 = openmc.Cell(fill=mat_clad, region=+surf_fuel & -surf_clad)
        c_he_hot2 = openmc.Cell(fill=mat_cool, region=+surf_clad & -surf_hot)
        c_part2 = openmc.Cell(fill=mat_clad, region=+surf_hot & -surf_part)
        c_he_cold2 = openmc.Cell(fill=mat_cool, region=+surf_part & -surf_cold)
        c_al2 = openmc.Cell(fill=mat_al2o3, region=+surf_cold & -surf_liner_da)
        c_mod2 = openmc.Cell(fill=mat_mod, region=+surf_liner_da)
        u_pin_out.add_cells([c_f2, c_cl2, c_he_hot2, c_part2, c_he_cold2, c_al2, c_mod2])

    # Hexagonal Lattice definition with orientation 'y'
    v1 = np.array([pitch, 0.0])
    v2 = np.array([pitch / 2.0, pitch * np.sqrt(3.0) / 2.0])

    ring_universes = []
    for r in range(n_rings, -1, -1):
        is_inner = (r <= 2)
        u_fill = u_pin_in if is_inner else u_pin_out
        ring_u = []
        if r == 0:
            ring_u.append(u_fill)
        else:
            corners = [r*v2, r*(v2-v1), -r*v1, -r*v2, -r*(v2-v1), r*v1]
            for corn in range(6):
                p_start = corners[corn]
                p_end = corners[(corn + 1) % 6]
                for step in range(r):
                    ring_u.append(u_fill)
        ring_universes.append(ring_u)

    hex_lat = openmc.HexLattice(name="Core_Hex_Lattice")
    hex_lat.center = (0.0, 0.0)
    hex_lat.pitch = (pitch,)
    hex_lat.orientation = 'y'
    hex_lat.universes = ring_universes
    # Outer background of lattice: YHx moderator
    u_mod_bg = openmc.Universe(name="Mod_Background")
    u_mod_bg.add_cell(openmc.Cell(fill=mat_mod))
    hex_lat.outer = u_mod_bg

    # Core active cell
    cell_core = openmc.Cell(name="Active_Core", fill=hex_lat,
                            region=-surf_cyl_core & +surf_z_bot_core & -surf_z_top_core)

    # -------------------------------------------------------------------------
    # 3. CONTROL DRUMS & RADIAL REFLECTOR
    # -------------------------------------------------------------------------
    n_drums = spec['drums']
    t_b4c = spec['t_b4c_cm']
    arc_b4c = spec['arc_b4c_deg']

    # Drum universe defined at origin (0, 0)
    cyl_drum_outer = openmc.ZCylinder(r=r_drum)
    cyl_drum_inner = openmc.ZCylinder(r=r_drum - t_b4c)
    half_sec = math.radians(arc_b4c / 2.0)
    plane_sec1 = openmc.Plane(a=math.sin(half_sec), b=-math.cos(half_sec), c=0.0, d=0.0)
    plane_sec2 = openmc.Plane(a=math.sin(half_sec), b=+math.cos(half_sec), c=0.0, d=0.0)
    plane_front = openmc.XPlane(x0=0.0)
    abs_reg = -cyl_drum_outer & +cyl_drum_inner & +plane_sec1 & +plane_sec2 & +plane_front
    refl_drum_reg = -cyl_drum_outer & ~abs_reg
    u_drum = openmc.Universe(cells=[
        openmc.Cell(fill=mat_b4c, region=abs_reg),
        openmc.Cell(fill=mat_refl, region=refl_drum_reg)
    ])

    drum_cells = []
    drum_union = None
    for k in range(n_drums):
        phi_k = k * 2.0 * math.pi / n_drums
        xk = R_drum_center * math.cos(phi_k)
        yk = R_drum_center * math.sin(phi_k)
        cyl_k = openmc.ZCylinder(x0=xk, y0=yk, r=r_drum)
        rot_deg = math.degrees(phi_k) + theta_deg
        c_dk = openmc.Cell(fill=u_drum, region=-cyl_k & +surf_z_bot_core & -surf_z_top_core)
        c_dk.translation = (xk, yk, 0.0)
        c_dk.rotation = (0.0, 0.0, rot_deg)
        drum_cells.append(c_dk)
        drum_union = -cyl_k if drum_union is None else (drum_union | -cyl_k)

    # Radial Reflector Cell (excluding drum cylinders)
    refl_rad_region = +surf_cyl_core & -surf_cyl_ref & +surf_z_bot_core & -surf_z_top_core & ~drum_union
    cell_refl_rad = openmc.Cell(name="Radial_Reflector", fill=mat_refl, region=refl_rad_region)

    # Axial Reflectors
    cell_refl_bot = openmc.Cell(name="Bottom_Axial_Reflector", fill=mat_refl,
                                region=-surf_cyl_ref & +surf_z_bot_tot & -surf_z_bot_core)
    cell_refl_top = openmc.Cell(name="Top_Axial_Reflector", fill=mat_refl,
                                region=-surf_cyl_ref & +surf_z_top_core & -surf_z_top_tot)

    root_universe = openmc.Universe(cells=[cell_core, cell_refl_rad, cell_refl_bot, cell_refl_top] + drum_cells)
    geometry = openmc.Geometry(root_universe)

    # -------------------------------------------------------------------------
    # 4. SETTINGS
    # -------------------------------------------------------------------------
    settings = openmc.Settings()
    settings.batches = batches
    settings.inactive = inactive
    settings.particles = particles
    settings.seed = seed
    settings.temperature = {'method': 'interpolation', 'tolerance': 10.0}

    # Initial spatial source in fuel region
    source_box = openmc.stats.Box(
        [-R_core * 0.7, -R_core * 0.7, -H_core * 0.4],
        [R_core * 0.7, R_core * 0.7, H_core * 0.4]
    )
    settings.source = openmc.IndependentSource(space=source_box)

    model = openmc.Model(materials=materials, geometry=geometry, settings=settings)
    return model

if __name__ == "__main__":
    from catalog import get_candidates
    c = get_candidates()[0]
    print(f"Testing model generation for {c['id']}...")
    model = build_model(c, theta_deg=0.0, particles=2000, batches=20, inactive=5)
    model.export_to_xml()
    print("Model XML files exported successfully!")
