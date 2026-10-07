#!/usr/bin/env python3
"""
build_reference_model.py

Frozen reference model generator for active reference reactor:
V7_CAND_BALANCED_ROBUST (Opt-17 Hexagonal Array, 217 pins, 400 kWth).
"""

import os
import math
import numpy as np
import openmc

# Global paths
OPENMC_EXEC = "/home/aleri/miniforge3/envs/openmc-env/bin/openmc"
CROSS_SECTIONS_PATH = "/home/aleri/wmp_data/cross_sections_wmp.xml"
CHAIN_PATH = "/home/aleri/OpenMC_projects/first_test/data/chain_nep.xml"

# Material densities [g/cm3]
RHO_UN = 13.50
RHO_CLAD = 8.57
RHO_YH2 = 4.30
RHO_BEO = 2.90
RHO_B4C = 2.52

# Thermophysical constants
P_NOM = 2.50e6          # 2.5 MPa
R_GAS = 8314.462618
M_HE = 4.002602
N_PINS = 217
DRUM_COUNT = 10


def rho_helium(P, T):
    """Density of pure helium [kg/m3]."""
    return (P * M_HE) / (R_GAS * max(T, 100.0))


def build_reference_openmc_model(run_dir,
                                 r_fuel=0.580,
                                 clad_th=0.048,
                                 chan_gap=0.084,
                                 pitch=2.270,
                                 core_h=42.45,
                                 refl_ax_th=9.275,
                                 refl_rad_th=9.80,
                                 enr_inner=19.30,
                                 enr_outer=19.75,
                                 drum_r=3.95,
                                 abs_th=1.50,
                                 abs_angle=148.0,
                                 drum_clearance=0.180,
                                 drum_angles=None,
                                 T_fuel_grid=None,
                                 T_clad_grid=None,
                                 T_cool_grid=None,
                                 T_mod=965.0,
                                 T_refl=900.0,
                                 T_fuel_offset=0.0,
                                 T_mod_offset=0.0,
                                 cool_density_mult=1.0,
                                 n_axial_slices=10,
                                 particles=50000,
                                 batches=35,
                                 inactive=10,
                                 seed=42,
                                 tallies_enabled=False):
    """Exports frozen OpenMC input XML files for V7_CAND_BALANCED_ROBUST."""
    os.makedirs(run_dir, exist_ok=True)
    openmc.reset_auto_ids()
    openmc.config['cross_sections'] = CROSS_SECTIONS_PATH

    if drum_angles is None:
        drum_angles = [0.0] * DRUM_COUNT

    if T_fuel_grid is None:
        T_fuel_grid = np.full((N_PINS, n_axial_slices), 1235.0)
    if T_clad_grid is None:
        T_clad_grid = np.full((N_PINS, n_axial_slices), 1220.0)
    if T_cool_grid is None:
        T_cool_grid = np.full((N_PINS, n_axial_slices), 1050.0)

    r_c = r_fuel + clad_th
    r_ch = r_c + chan_gap
    r_add = 0.150
    r_add_pos = 0.910 * (pitch / 2.270)

    core_r = 19.45 * (pitch / 2.270)
    refl_rad_outer_r = core_r + refl_rad_th
    total_h = core_h + 2.0 * refl_ax_th
    ax_duct_r = 0.860
    drum_rad_pos = core_r + drum_clearance + drum_r

    T_fuel_grid_eff = T_fuel_grid + T_fuel_offset
    T_mod_eff = T_mod + T_mod_offset

    # 1. Materials
    # Fuel Inner (Rings 0-2)
    f_in = openmc.Material(name=f"UN Fuel Inner ({enr_inner:.2f}%)")
    f_in.set_density("g/cm3", RHO_UN)
    u235_at_in = (enr_inner / 235.0) / (enr_inner / 235.0 + (100.0 - enr_inner) / 238.0)
    f_in.add_nuclide("U235", 0.5 * u235_at_in, "ao")
    f_in.add_nuclide("U238", 0.5 * (1.0 - u235_at_in), "ao")
    f_in.add_nuclide("N14", 0.5 * 0.9963, "ao")
    f_in.add_nuclide("N15", 0.5 * 0.0037, "ao")
    f_in.temperature = float(np.mean(T_fuel_grid_eff[:19, :]))

    # Fuel Outer (Rings 3-8)
    f_out = openmc.Material(name=f"UN Fuel Outer ({enr_outer:.2f}%)")
    f_out.set_density("g/cm3", RHO_UN)
    u235_at_out = (enr_outer / 235.0) / (enr_outer / 235.0 + (100.0 - enr_outer) / 238.0)
    f_out.add_nuclide("U235", 0.5 * u235_at_out, "ao")
    f_out.add_nuclide("U238", 0.5 * (1.0 - u235_at_out), "ao")
    f_out.add_nuclide("N14", 0.5 * 0.9963, "ao")
    f_out.add_nuclide("N15", 0.5 * 0.0037, "ao")
    f_out.temperature = float(np.mean(T_fuel_grid_eff[19:, :]))

    # Cladding Nb-1Zr
    clad = openmc.Material(name="Nb-1Zr Clad")
    clad.set_density("g/cm3", RHO_CLAD)
    clad.add_element("Nb", 0.99, "wo")
    clad.add_element("Zr", 0.01, "wo")
    clad.temperature = float(np.mean(T_clad_grid))

    # Helium Coolant slices (10 axial slices)
    cool_axial_mats = []
    for s in range(n_axial_slices):
        T_cool_s = float(np.mean(T_cool_grid[:, s]))
        rho_he_si = (rho_helium(P_NOM, T_cool_s) / 1000.0) * cool_density_mult
        m_cool_s = openmc.Material(name=f"He Coolant S{s}")
        m_cool_s.set_density("g/cm3", max(rho_he_si, 1e-7))
        m_cool_s.add_nuclide("He4", 1.0, "ao")
        m_cool_s.temperature = T_cool_s
        cool_axial_mats.append(m_cool_s)

    # Plenum Coolant
    rho_plenum = (rho_helium(P_NOM, 1050.0) / 1000.0) * cool_density_mult
    cool_plenum = openmc.Material(name="He Coolant Plenum")
    cool_plenum.set_density("g/cm3", max(rho_plenum, 1e-7))
    cool_plenum.add_nuclide("He4", 1.0, "ao")
    cool_plenum.temperature = 1050.0

    # Moderator YH2.0 (stoichiometric canned hydride)
    mod = openmc.Material(name="Canned YH2.0")
    mod.set_density("g/cm3", RHO_YH2)
    mod.add_nuclide("Y89", 1.0, "ao")
    mod.add_nuclide("H1", 2.0, "ao")
    mod.add_s_alpha_beta("c_H_in_ZrH")
    mod.temperature = T_mod_eff

    # Reflector BeO
    refl = openmc.Material(name="BeO Reflector")
    refl.set_density("g/cm3", RHO_BEO)
    refl.add_nuclide("Be9", 1.0, "ao")
    refl.add_nuclide("O16", 1.0, "ao")
    refl.add_s_alpha_beta("c_Be_in_BeO")
    refl.add_s_alpha_beta("c_O_in_BeO")
    refl.temperature = T_refl

    # Absorber B4C (90% enriched B-10)
    absorb = openmc.Material(name="B4C Absorber")
    absorb.set_density("g/cm3", RHO_B4C)
    absorb.add_nuclide("B10", 4.0 * 0.90, "ao")
    absorb.add_nuclide("B11", 4.0 * 0.10, "ao")
    absorb.add_nuclide("C0", 1.0, "ao")
    absorb.temperature = T_refl

    materials = openmc.Materials([f_in, f_out, clad, cool_plenum, mod, refl, absorb] + cool_axial_mats)

    # 2. Geometry
    z_cb = openmc.ZPlane(z0=-core_h / 2.0)
    z_ct = openmc.ZPlane(z0=+core_h / 2.0)
    z_rb = openmc.ZPlane(z0=-total_h / 2.0, boundary_type="vacuum")
    z_rt = openmc.ZPlane(z0=+total_h / 2.0, boundary_type="vacuum")

    z_slice_edges = np.linspace(-core_h / 2.0, +core_h / 2.0, n_axial_slices + 1)
    z_slice_planes = [openmc.ZPlane(z0=z_val) for z_val in z_slice_edges]

    cyl_f = openmc.ZCylinder(r=r_fuel)
    cyl_c = openmc.ZCylinder(r=r_c)
    cyl_ch = openmc.ZCylinder(r=r_ch)

    aux_cyls = []
    for k in range(6):
        ang = math.radians(k * 60.0)
        aux_cyls.append(openmc.ZCylinder(x0=r_add_pos * math.cos(ang),
                                         y0=r_add_pos * math.sin(ang),
                                         r=r_add))

    cyl_ax = openmc.ZCylinder(r=ax_duct_r)

    fuel_cells_map = {}
    cool_cells_map = {}
    clad_cells_map = {}

    def make_pin_unique(pin_idx, f_mat):
        pin_cells = []
        for s_idx in range(n_axial_slices):
            z_low = z_slice_planes[s_idx]
            z_high = z_slice_planes[s_idx + 1]
            T_f_cell = float(T_fuel_grid_eff[pin_idx, s_idx])
            T_c_cell = float(T_clad_grid[pin_idx, s_idx])
            T_cool_cell = float(T_cool_grid[pin_idx, s_idx])
            cool_mat = cool_axial_mats[s_idx]

            c_f = openmc.Cell(name=f"Fuel_P{pin_idx}_S{s_idx}", fill=f_mat, region=-cyl_f & +z_low & -z_high)
            c_f.temperature = T_f_cell
            pin_cells.append(c_f)
            fuel_cells_map[(pin_idx, s_idx)] = c_f

            c_c = openmc.Cell(name=f"Clad_P{pin_idx}_S{s_idx}", fill=clad, region=+cyl_f & -cyl_c & +z_low & -z_high)
            c_c.temperature = T_c_cell
            pin_cells.append(c_c)
            clad_cells_map[(pin_idx, s_idx)] = c_c

            c_ann = openmc.Cell(name=f"CoolAnn_P{pin_idx}_S{s_idx}", fill=cool_mat, region=+cyl_c & -cyl_ch & +z_low & -z_high)
            c_ann.temperature = T_cool_cell
            pin_cells.append(c_ann)
            cool_cells_map[(pin_idx, s_idx)] = c_ann

        aux_cells = []
        aux_union = None
        for c in aux_cyls:
            aux_cells.append(openmc.Cell(fill=cool_plenum, region=-c & +z_cb & -z_ct))
            aux_union = -c if aux_union is None else (aux_union | -c)

        mod_reg = +cyl_ch & +z_cb & -z_ct
        if aux_union is not None:
            mod_reg = mod_reg & ~aux_union
        c_m = openmc.Cell(fill=mod, region=mod_reg)

        # Axial reflector extensions
        c_db = openmc.Cell(fill=cool_plenum, region=-cyl_ax & +z_rb & -z_cb)
        c_dt = openmc.Cell(fill=cool_plenum, region=-cyl_ax & +z_ct & -z_rt)

        aux_db = []
        aux_dt = []
        aux_ax_union = None
        for c in aux_cyls:
            aux_db.append(openmc.Cell(fill=cool_plenum, region=-c & +z_rb & -z_cb))
            aux_dt.append(openmc.Cell(fill=cool_plenum, region=-c & +z_ct & -z_rt))
            aux_ax_union = -c if aux_ax_union is None else (aux_ax_union | -c)

        c_bb = openmc.Cell(fill=refl, region=+cyl_ax & +z_rb & -z_cb & ~aux_ax_union)
        c_bb.temperature = T_refl
        c_bt = openmc.Cell(fill=refl, region=+cyl_ax & +z_ct & -z_rt & ~aux_ax_union)
        c_bt.temperature = T_refl

        u = openmc.Universe(cells=pin_cells + [c_m, c_db, c_bb, c_dt, c_bt] + aux_cells + aux_db + aux_dt)
        return u

    # Hexagonal Lattice
    v1 = np.array([pitch, 0.0])
    v2 = np.array([pitch / 2.0, pitch * np.sqrt(3.0) / 2.0])

    pin_meta = []
    pin_idx_counter = 0

    rings = []
    for r in range(8, -1, -1):
        is_inner = (r <= 2)
        f_mat = f_in if is_inner else f_out
        ring_u = []
        if r == 0:
            u_p = make_pin_unique(pin_idx_counter, f_mat)
            ring_u.append(u_p)
            pin_meta.append((pin_idx_counter, 0, 0.0, 0.0, is_inner))
            pin_idx_counter += 1
        else:
            corners = [r*v2, r*(v2-v1), -r*v1, -r*v2, -r*(v2-v1), r*v1]
            for corn in range(6):
                p_start = corners[corn]
                p_end = corners[(corn + 1) % 6]
                for step in range(r):
                    pt = p_start + (p_end - p_start) * (step / r)
                    u_p = make_pin_unique(pin_idx_counter, f_mat)
                    ring_u.append(u_p)
                    pin_meta.append((pin_idx_counter, r, float(pt[0]), float(pt[1]), is_inner))
                    pin_idx_counter += 1
        rings.append(ring_u)

    u_mod_buf = openmc.Universe(cells=[
        openmc.Cell(fill=mod, region=+z_cb & -z_ct),
        openmc.Cell(fill=refl, region=+z_rb & -z_cb),
        openmc.Cell(fill=refl, region=+z_ct & -z_rt)
    ])

    hex_lat = openmc.HexLattice()
    hex_lat.center = (0.0, 0.0)
    hex_lat.pitch = [pitch]
    hex_lat.orientation = "y"
    hex_lat.universes = rings
    hex_lat.outer = u_mod_buf

    cyl_core = openmc.ZCylinder(r=core_r)
    cyl_refl_outer = openmc.ZCylinder(r=refl_rad_outer_r, boundary_type="vacuum")

    # 10 Control Drums
    cyl_drum_outer = openmc.ZCylinder(r=drum_r)
    cyl_drum_inner = openmc.ZCylinder(r=drum_r - abs_th)
    half_sec = np.radians(abs_angle / 2.0)
    plane_sec1 = openmc.Plane(a=np.sin(half_sec), b=-np.cos(half_sec), c=0.0, d=0.0)
    plane_sec2 = openmc.Plane(a=np.sin(half_sec), b=+np.cos(half_sec), c=0.0, d=0.0)
    plane_front = openmc.XPlane(x0=0.0)
    abs_reg = -cyl_drum_outer & +cyl_drum_inner & +plane_sec1 & +plane_sec2 & +plane_front
    beo_drum_reg = -cyl_drum_outer & ~abs_reg
    u_drum = openmc.Universe(cells=[
        openmc.Cell(fill=absorb, region=abs_reg),
        openmc.Cell(fill=refl, region=beo_drum_reg)
    ])

    drum_cells = []
    drum_regs = []
    for k in range(DRUM_COUNT):
        phi_k = k * 2.0 * np.pi / DRUM_COUNT
        xk = drum_rad_pos * math.cos(phi_k)
        yk = drum_rad_pos * math.sin(phi_k)
        cyl_k = openmc.ZCylinder(x0=xk, y0=yk, r=drum_r)
        rot_deg = k * (360.0 / DRUM_COUNT) + drum_angles[k]
        c_dk = openmc.Cell(fill=u_drum, region=-cyl_k & +z_cb & -z_ct)
        c_dk.translation = (xk, yk, 0.0)
        c_dk.rotation = (0.0, 0.0, rot_deg)
        drum_cells.append(c_dk)
        drum_regs.append(-cyl_k)

    drum_union = drum_regs[0]
    for reg in drum_regs[1:]:
        drum_union = drum_union | reg

    cell_core = openmc.Cell(fill=hex_lat, region=-cyl_core & +z_rb & -z_rt)
    cell_refl_rad = openmc.Cell(fill=refl, region=+cyl_core & -cyl_refl_outer & +z_cb & -z_ct & ~drum_union)
    cell_refl_bot = openmc.Cell(fill=refl, region=+cyl_core & -cyl_refl_outer & +z_rb & -z_cb)
    cell_refl_top = openmc.Cell(fill=refl, region=+cyl_core & -cyl_refl_outer & +z_ct & -z_rt)

    root = openmc.Universe(cells=[cell_core, cell_refl_rad, cell_refl_bot, cell_refl_top] + drum_cells)
    geometry = openmc.Geometry(root)

    # Tallies
    tallies = None
    if tallies_enabled:
        fuel_cell_ids = [c.id for c in fuel_cells_map.values()]
        cell_filter = openmc.CellFilter(fuel_cell_ids)
        t_pin_heat = openmc.Tally(name="fuel_cell_heating")
        t_pin_heat.filters = [cell_filter]
        t_pin_heat.scores = ["heating", "fission"]
        tallies = openmc.Tallies([t_pin_heat])

    # Settings
    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    settings.temperature = {"multipole": True, "method": "nearest", "tolerance": 1500.0}
    settings.seed = seed

    entropy_mesh = openmc.RegularMesh()
    entropy_mesh.lower_left = (-core_r, -core_r, -core_h / 2.0)
    entropy_mesh.upper_right = (core_r, core_r, core_h / 2.0)
    entropy_mesh.dimension = (10, 10, min(n_axial_slices, 10))
    settings.entropy_mesh = entropy_mesh

    space = openmc.stats.Box((-core_r * 0.6, -core_r * 0.6, -core_h * 0.3),
                            (core_r * 0.6, core_r * 0.6, core_h * 0.3))
    settings.source = openmc.IndependentSource(space=space)
    settings.sourcepoint = {"write": True}

    materials.export_to_xml(os.path.join(run_dir, "materials.xml"))
    geometry.export_to_xml(os.path.join(run_dir, "geometry.xml"))
    settings.export_to_xml(os.path.join(run_dir, "settings.xml"))
    if tallies is not None:
        tallies.export_to_xml(os.path.join(run_dir, "tallies.xml"))

    return {
        "materials": materials,
        "geometry": geometry,
        "settings": settings,
        "tallies": tallies,
        "pin_meta": pin_meta,
        "fuel_cells_map": fuel_cells_map,
        "z_slice_edges": z_slice_edges,
        "n_axial_slices": n_axial_slices,
        "core_r": core_r,
        "refl_rad_outer_r": refl_rad_outer_r,
        "total_h": total_h
    }
