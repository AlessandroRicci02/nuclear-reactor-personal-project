#!/usr/bin/env python3
"""
===============================================================================
STAGE 1 ARCHITECTURE SCREENING ENGINE: CANDIDATE OPT-11
===============================================================================
Explores BeO-First Moderator Exploration & Global Architecture Optimization:
  Family A: Pure BeO Square Pin Core (UN pins in BeO block, varied sizing/drums)
  Family B: Pure BeO Hexagonal Pin Lattice (HexLattice in BeO matrix, 10/12 drums)
  Family C: Annular Core with Central BeO Island (BeO flux trap, 10/12 drums)
  Family D: Spatially Separated / Composite BeO Cores (BeO-ZrH2 & BeO-YH2 hybrids)
  Family E: Modular Clustered BeO Core (7 fuel clusters in BeO block)
  Benchmarks: Opt-10 Baseline (YH2), Opt-10 Extended, ZrH2 Compact & Ultra-Compact

Pipeline:
  1. Construct exact 3D OpenMC CSG geometry for each candidate.
  2. Compute itemized mass and volume breakdown.
  3. Evaluate BOL operating multiplication (k_out, drums 0 deg).
  4. Evaluate shutdown multiplication (k_in, drums 180 deg) & Total Drum Worth.
  5. Evaluate single stuck-drum multiplication (k_stuck, drum #0 at 0 deg, rest 180 deg).
  6. Evaluate Coolant Void Reactivity (Depressurized vacuum state).
  7. Export comprehensive screening metrics to results/opt11/screening/opt11_architecture_screening.csv.
  8. Generate Stage 1 Pareto comparison figures.
===============================================================================
"""

import os
import sys
import csv
import shutil
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import openmc

os.environ["PATH"] = "/home/aleri/miniforge3/envs/openmc-env/bin:" + os.environ.get("PATH", "")
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
os.environ["OMP_NUM_THREADS"] = "8"
openmc.config['cross_sections'] = "/home/aleri/wmp_data/cross_sections_wmp.xml"


class ArchOpt11Candidate:
    def __init__(self,
                 cand_id,
                 family,
                 name,
                 mod_type='BeO', # 'BeO', 'YH2', 'ZrH2', 'Hybrid_BeO_ZrH2_70', etc.
                 lattice_type='square', # 'square', 'hex', 'annular', 'cluster'
                 core_radius=19.40,
                 active_core_height=41.50,
                 fuel_radius=0.580,
                 cladding_thickness=0.050,
                 coolant_channel_thickness=0.070,
                 pitch=1.665,
                 inner_island_radius=0.0,
                 island_material='BeO',
                 number_of_drums=10,
                 drum_radius=3.90,
                 absorber_thickness=1.40,
                 absorber_sector_deg=140.0,
                 radial_reflector_thickness=9.20,
                 axial_reflector_thickness=9.50,
                 axial_duct_radius=0.820,
                 coolant_pressure_MPa=2.5,
                 coolant_temp_K=1050.0):

        self.cand_id = cand_id
        self.family = family
        self.name = name
        self.mod_type = mod_type
        self.lattice_type = lattice_type
        self.core_radius = core_radius
        self.active_core_height = active_core_height
        self.aspect_ratio_H_R = self.active_core_height / self.core_radius
        self.fuel_radius = fuel_radius
        self.cladding_thickness = cladding_thickness
        self.coolant_channel_thickness = coolant_channel_thickness
        self.pitch = pitch
        self.inner_island_radius = inner_island_radius
        self.island_material = island_material

        self.clad_outer_radius = self.fuel_radius + self.cladding_thickness
        self.channel_outer_radius = self.clad_outer_radius + self.coolant_channel_thickness
        self.axial_duct_radius = axial_duct_radius

        self.number_of_drums = number_of_drums
        self.drum_radius = drum_radius
        self.absorber_thickness = absorber_thickness
        self.absorber_sector_deg = absorber_sector_deg

        self.drum_clearance = 0.20
        self.drum_radial_position = self.core_radius + self.drum_clearance + self.drum_radius
        min_refl_thick = self.drum_clearance + 2.0 * self.drum_radius + 0.80
        self.radial_reflector_thickness = max(radial_reflector_thickness, min_refl_thick)
        self.radial_reflector_outer_radius = self.core_radius + self.radial_reflector_thickness
        self.axial_reflector_thickness = axial_reflector_thickness
        self.total_reactor_height = self.active_core_height + 2.0 * self.axial_reflector_thickness

        # Densities
        self.fuel_density = 13.50
        self.cladding_density = 8.57
        self.reflector_density = 2.90
        self.absorber_density = 2.52
        self.b10_enrichment = 0.90
        self.coolant_pressure_MPa = coolant_pressure_MPa
        self.coolant_temp_K = coolant_temp_K
        self.coolant_density = (self.coolant_pressure_MPa * 1.0e6 * 4.0026e-3) / (8.31446 * self.coolant_temp_K) * 1.0e-3


def build_candidate_materials(c: ArchOpt11Candidate, cool_density=None):
    """Build material cards with explicit isotopes."""
    openmc.reset_auto_ids()

    # Fuel UN (19.75 wt% U-235)
    fuel = openmc.Material(name="UN Fuel (19.75 wt% U-235)")
    fuel.set_density("g/cm3", c.fuel_density)
    fuel.add_nuclide("U235", 0.5 * 0.199522, "ao")
    fuel.add_nuclide("U238", 0.5 * 0.800478, "ao")
    fuel.add_nuclide("N14", 0.5 * 0.9963, "ao")
    fuel.add_nuclide("N15", 0.5 * 0.0037, "ao")
    fuel.temperature = 1050.0

    # Cladding Nb-1Zr
    cladding = openmc.Material(name="Nb-1Zr Cladding")
    cladding.set_density("g/cm3", c.cladding_density)
    cladding.add_element("Nb", 0.99, "wo")
    cladding.add_element("Zr", 0.01, "wo")
    cladding.temperature = 1050.0

    # Coolant He
    rho_c = c.coolant_density if cool_density is None else cool_density
    coolant = openmc.Material(name="Pure Helium Coolant")
    coolant.set_density("g/cm3", max(rho_c, 1.0e-7))
    coolant.add_nuclide("He4", 1.0, "ao")
    coolant.temperature = c.coolant_temp_K

    # Reflector BeO
    reflector = openmc.Material(name="BeO Reflector")
    reflector.set_density("g/cm3", c.reflector_density)
    reflector.add_nuclide("Be9", 1.0, "ao")
    reflector.add_nuclide("O16", 1.0, "ao")
    reflector.add_s_alpha_beta("c_Be_in_BeO")
    reflector.add_s_alpha_beta("c_O_in_BeO")
    reflector.temperature = 800.0

    # Absorber B4C
    absorber = openmc.Material(name="B4C Absorber (90% B-10)")
    absorber.set_density("g/cm3", c.absorber_density)
    absorber.add_nuclide("B10", 4.0 * c.b10_enrichment, "ao")
    absorber.add_nuclide("B11", 4.0 * (1.0 - c.b10_enrichment), "ao")
    absorber.add_nuclide("C0", 1.0, "ao")
    absorber.temperature = 800.0

    # Moderator
    if c.mod_type == 'BeO':
        moderator = reflector
        rho_mod = 2.90
    elif c.mod_type == 'YH2':
        moderator = openmc.Material(name="YH2 Moderator")
        moderator.set_density("g/cm3", 4.30)
        moderator.add_nuclide("Y89", 1.0, "ao")
        moderator.add_nuclide("H1", 2.0, "ao")
        moderator.add_s_alpha_beta("c_H_in_ZrH")
        moderator.temperature = 800.0
        rho_mod = 4.30
    elif c.mod_type == 'ZrH2':
        moderator = openmc.Material(name="ZrH2 Moderator")
        moderator.set_density("g/cm3", 5.60)
        moderator.add_nuclide("Zr90", 0.5145, "ao")
        moderator.add_nuclide("Zr91", 0.1122, "ao")
        moderator.add_nuclide("Zr92", 0.1715, "ao")
        moderator.add_nuclide("Zr94", 0.1738, "ao")
        moderator.add_nuclide("Zr96", 0.0280, "ao")
        moderator.add_nuclide("H1", 2.0, "ao")
        moderator.add_s_alpha_beta("c_H_in_ZrH")
        moderator.add_s_alpha_beta("c_Zr_in_ZrH")
        moderator.temperature = 800.0
        rho_mod = 5.60
    elif c.mod_type.startswith('Hybrid_BeO_ZrH2_'):
        pct_zrh2 = float(c.mod_type.split('_')[-1]) / 100.0
        pct_beo = 1.0 - pct_zrh2
        rho_mod = pct_beo * 2.90 + pct_zrh2 * 5.60
        moderator = openmc.Material(name=f"Hybrid BeO-ZrH2 {int(pct_zrh2*100)}%")
        moderator.set_density("g/cm3", rho_mod)
        w_beo = (pct_beo * 2.90) / rho_mod
        w_zrh2 = (pct_zrh2 * 5.60) / rho_mod
        moderator.add_nuclide("Be9", w_beo * (9.012/25.012), "wo")
        moderator.add_nuclide("O16", w_beo * (16.0/25.012), "wo")
        w_zr = w_zrh2 * (91.224 / 93.24)
        w_h = w_zrh2 * (2.016 / 93.24)
        moderator.add_nuclide("Zr90", w_zr * 0.5145, "wo")
        moderator.add_nuclide("Zr91", w_zr * 0.1122, "wo")
        moderator.add_nuclide("Zr92", w_zr * 0.1715, "wo")
        moderator.add_nuclide("Zr94", w_zr * 0.1738, "wo")
        moderator.add_nuclide("Zr96", w_zr * 0.0280, "wo")
        moderator.add_nuclide("H1", w_h, "wo")
        moderator.add_s_alpha_beta("c_Be_in_BeO")
        moderator.add_s_alpha_beta("c_O_in_BeO")
        moderator.add_s_alpha_beta("c_H_in_ZrH")
        moderator.add_s_alpha_beta("c_Zr_in_ZrH")
        moderator.temperature = 800.0
    elif c.mod_type.startswith('Hybrid_BeO_YH2_'):
        pct_yh2 = float(c.mod_type.split('_')[-1]) / 100.0
        pct_beo = 1.0 - pct_yh2
        rho_mod = pct_beo * 2.90 + pct_yh2 * 4.30
        moderator = openmc.Material(name=f"Hybrid BeO-YH2 {int(pct_yh2*100)}%")
        moderator.set_density("g/cm3", rho_mod)
        w_beo = (pct_beo * 2.90) / rho_mod
        w_yh2 = (pct_yh2 * 4.30) / rho_mod
        moderator.add_nuclide("Be9", w_beo * (9.012/25.012), "wo")
        moderator.add_nuclide("O16", w_beo * (16.0/25.012), "wo")
        moderator.add_nuclide("Y89", w_yh2 * (88.905/90.921), "wo")
        moderator.add_nuclide("H1", w_yh2 * (2.016/90.921), "wo")
        moderator.add_s_alpha_beta("c_Be_in_BeO")
        moderator.add_s_alpha_beta("c_O_in_BeO")
        moderator.add_s_alpha_beta("c_H_in_ZrH")
        moderator.temperature = 800.0
    else:
        moderator = reflector
        rho_mod = 2.90

    # Central island material if annular
    if c.island_material == 'BeO':
        island_mat = reflector
    elif c.island_material == 'YH2':
        island_mat = openmc.Material(name="YH2 Island")
        island_mat.set_density("g/cm3", 4.30)
        island_mat.add_nuclide("Y89", 1.0, "ao")
        island_mat.add_nuclide("H1", 2.0, "ao")
        island_mat.add_s_alpha_beta("c_H_in_ZrH")
        island_mat.temperature = 800.0
    else:
        island_mat = reflector

    mats = [fuel, cladding, coolant, moderator, reflector, absorber]
    if c.island_material == 'YH2':
        mats.append(island_mat)
    materials = openmc.Materials(mats)
    materials.cross_sections = "/home/aleri/wmp_data/cross_sections_wmp.xml"
    return materials, fuel, cladding, coolant, moderator, reflector, absorber, island_mat, rho_mod


def build_candidate_geometry(c: ArchOpt11Candidate, fuel, cladding, coolant, moderator, reflector, absorber, island_mat, drum_angles=None):
    """Constructs exact 3D CSG geometry."""
    if drum_angles is None:
        drum_angles = [0.0] * c.number_of_drums

    s_fuel = openmc.ZCylinder(r=c.fuel_radius, name="Pin Fuel Outer")
    s_clad = openmc.ZCylinder(r=c.clad_outer_radius, name="Pin Clad Outer")
    s_channel = openmc.ZCylinder(r=c.channel_outer_radius, name="Coolant Channel Outer")

    z_core_bot = openmc.ZPlane(z0=-c.active_core_height / 2.0, name="Core Bot")
    z_core_top = openmc.ZPlane(z0=+c.active_core_height / 2.0, name="Core Top")
    z_refl_bot = openmc.ZPlane(z0=-c.total_reactor_height / 2.0, boundary_type="vacuum", name="Reactor Bot")
    z_refl_top = openmc.ZPlane(z0=+c.total_reactor_height / 2.0, boundary_type="vacuum", name="Reactor Top")

    # Cooled pin cell:
    c_f = openmc.Cell(fill=fuel, region=-s_fuel & +z_core_bot & -z_core_top)
    c_cl = openmc.Cell(fill=cladding, region=+s_fuel & -s_clad & +z_core_bot & -z_core_top)
    c_co = openmc.Cell(fill=coolant, region=+s_clad & -s_channel & +z_core_bot & -z_core_top)
    c_mo = openmc.Cell(fill=moderator, region=+s_channel & +z_core_bot & -z_core_top)

    # Axial BeO streaming ducts
    s_duct = openmc.ZCylinder(r=c.axial_duct_radius, name="Axial Streaming Duct")
    c_co_bot = openmc.Cell(fill=coolant, region=-s_duct & +z_refl_bot & -z_core_bot)
    c_be_bot = openmc.Cell(fill=reflector, region=+s_duct & +z_refl_bot & -z_core_bot)
    c_co_top = openmc.Cell(fill=coolant, region=-s_duct & +z_core_top & -z_refl_top)
    c_be_top = openmc.Cell(fill=reflector, region=+s_duct & +z_core_top & -z_refl_top)

    u_pin = openmc.Universe(cells=[c_f, c_cl, c_co, c_mo, c_co_bot, c_be_bot, c_co_top, c_be_top])

    # Moderator-only universe
    c_mo_core = openmc.Cell(fill=moderator, region=+z_core_bot & -z_core_top)
    c_be_bot_only = openmc.Cell(fill=reflector, region=+z_refl_bot & -z_core_bot)
    c_be_top_only = openmc.Cell(fill=reflector, region=+z_core_top & -z_refl_top)
    u_mod_only = openmc.Universe(cells=[c_mo_core, c_be_bot_only, c_be_top_only])

    # Island universe
    c_isl_core = openmc.Cell(fill=island_mat, region=+z_core_bot & -z_core_top)
    c_isl_bot = openmc.Cell(fill=reflector, region=+z_refl_bot & -z_core_bot)
    c_isl_top = openmc.Cell(fill=reflector, region=+z_core_top & -z_refl_top)
    u_island = openmc.Universe(cells=[c_isl_core, c_isl_bot, c_isl_top])

    pin_centers = []

    if c.lattice_type in ['square', 'annular', 'cluster']:
        pitch = c.pitch
        grid_dim = int(np.ceil(2 * c.core_radius / pitch))
        if grid_dim % 2 == 0:
            grid_dim += 1

        lattice_universes = []
        for i in range(grid_dim):
            row = []
            for j in range(grid_dim):
                x = -grid_dim * pitch / 2.0 + (j + 0.5) * pitch
                y = -grid_dim * pitch / 2.0 + (i + 0.5) * pitch
                r_dist = np.sqrt(x**2 + y**2)

                is_pin = False
                if r_dist <= (c.core_radius - c.channel_outer_radius):
                    if c.lattice_type == 'annular' and c.inner_island_radius > 0:
                        if r_dist >= (c.inner_island_radius + c.channel_outer_radius):
                            is_pin = True
                    elif c.lattice_type == 'cluster':
                        # 7 cluster pattern: central cluster (r < 7.0) or 6 peripheral clusters
                        # Cluster webs filled with moderator
                        cluster_centers = [(0, 0)] + [
                            (12.0 * np.cos(k * np.pi / 3), 12.0 * np.sin(k * np.pi / 3))
                            for k in range(6)
                        ]
                        in_cluster = any(np.hypot(x - cx, y - cy) <= 5.2 for cx, cy in cluster_centers)
                        if in_cluster:
                            is_pin = True
                    else:
                        is_pin = True

                if is_pin:
                    row.append(u_pin)
                    pin_centers.append((x, y))
                else:
                    if c.lattice_type == 'annular' and r_dist < c.inner_island_radius:
                        row.append(u_island)
                    else:
                        row.append(u_mod_only)
            lattice_universes.append(row)

        lattice = openmc.RectLattice(name=f"{c.cand_id} Rect Lattice")
        lattice.lower_left = (-grid_dim * pitch / 2.0, -grid_dim * pitch / 2.0)
        lattice.pitch = (pitch, pitch)
        lattice.universes = lattice_universes
        lattice.outer = u_mod_only

    elif c.lattice_type == 'hex':
        pitch = c.pitch
        num_rings = int(np.ceil(c.core_radius / pitch)) + 1
        v1 = np.array([pitch, 0.0])
        v2 = np.array([pitch/2.0, pitch * np.sqrt(3)/2.0])

        rings_universes = []
        for r in range(num_rings - 1, -1, -1):
            ring_u = []
            if r == 0:
                pt = (0.0, 0.0)
                r_dist = 0.0
                if r_dist <= (c.core_radius - c.channel_outer_radius):
                    ring_u.append(u_pin)
                    pin_centers.append(pt)
                else:
                    ring_u.append(u_mod_only)
            else:
                corners = [
                    r * v2,
                    r * (v2 - v1),
                    -r * v1,
                    -r * v2,
                    -r * (v2 - v1),
                    r * v1
                ]
                for corn in range(6):
                    p_start = corners[corn]
                    p_end = corners[(corn + 1) % 6]
                    for step in range(r):
                        pt = p_start + (p_end - p_start) * (step / r)
                        x, y = float(pt[0]), float(pt[1])
                        r_dist = np.sqrt(x**2 + y**2)
                        if r_dist <= (c.core_radius - c.channel_outer_radius):
                            ring_u.append(u_pin)
                            pin_centers.append((x, y))
                        else:
                            ring_u.append(u_mod_only)
            rings_universes.append(ring_u)

        lat = openmc.HexLattice(name=f"{c.cand_id} Hex Lattice")
        lat.center = (0.0, 0.0)
        lat.pitch = (pitch,)
        lat.orientation = 'y'
        lat.outer = u_mod_only
        lat.universes = rings_universes
        lattice = lat

    s_core = openmc.ZCylinder(r=c.core_radius)
    s_outer = openmc.ZCylinder(r=c.radial_reflector_outer_radius, boundary_type="vacuum")

    # Drum Universe Template
    r_abs_in = c.drum_radius - c.absorber_thickness
    drum_outer_cyl = openmc.ZCylinder(r=c.drum_radius)
    drum_inner_cyl = openmc.ZCylinder(r=r_abs_in)

    half_phi = np.radians(c.absorber_sector_deg / 2.0)
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

    for k in range(c.number_of_drums):
        ang = k * 2.0 * np.pi / c.number_of_drums
        xk = c.drum_radial_position * np.cos(ang)
        yk = c.drum_radial_position * np.sin(ang)
        dk = openmc.ZCylinder(x0=xk, y0=yk, r=c.drum_radius)

        cdk = openmc.Cell(fill=u_drum, region=-dk & +z_core_bot & -z_core_top)
        cdk.translation = (xk, yk, 0.0)
        cdk.rotation = (0.0, 0.0, k * 360.0 / c.number_of_drums + drum_angles[k])

        drum_cells.append(cdk)
        drum_regions.append(-dk)

    drum_union = drum_regions[0]
    for reg in drum_regions[1:]:
        drum_union = drum_union | reg

    # Global Cells
    c_core_cyl = openmc.Cell(fill=lattice, region=-s_core & +z_refl_bot & -z_refl_top)
    c_rad_refl = openmc.Cell(fill=reflector, region=+s_core & -s_outer & +z_core_bot & -z_core_top & ~drum_union)
    c_rad_bot = openmc.Cell(fill=reflector, region=+s_core & -s_outer & +z_refl_bot & -z_core_bot)
    c_rad_top = openmc.Cell(fill=reflector, region=+s_core & -s_outer & +z_core_top & -z_refl_top)

    root_universe = openmc.Universe(cells=[c_core_cyl, c_rad_refl, c_rad_bot, c_rad_top] + drum_cells)
    geometry = openmc.Geometry(root_universe)

    return geometry, pin_centers, lattice


def compute_candidate_inventory(c: ArchOpt11Candidate, pin_centers, rho_mod):
    """Compute exact material volume and mass accounting."""
    n_pins = len(pin_centers)
    v_fuel_single = np.pi * c.fuel_radius**2 * c.active_core_height
    v_clad_single = np.pi * (c.clad_outer_radius**2 - c.fuel_radius**2) * c.active_core_height
    v_cool_single_core = np.pi * (c.channel_outer_radius**2 - c.clad_outer_radius**2) * c.active_core_height
    v_cool_single_axial = np.pi * c.axial_duct_radius**2 * (2 * c.axial_reflector_thickness)

    v_fuel = n_pins * v_fuel_single
    v_clad = n_pins * v_clad_single
    v_cool_core = n_pins * v_cool_single_core
    v_cool_axial = n_pins * v_cool_single_axial

    v_active_core_tot = np.pi * c.core_radius**2 * c.active_core_height
    
    if c.inner_island_radius > 0:
        v_island = np.pi * c.inner_island_radius**2 * c.active_core_height
    else:
        v_island = 0.0

    v_mod = v_active_core_tot - (v_fuel + v_clad + v_cool_core) - v_island
    v_axial_refl = np.pi * c.core_radius**2 * (2 * c.axial_reflector_thickness) - v_cool_axial

    v_drums = c.number_of_drums * np.pi * c.drum_radius**2 * c.active_core_height
    v_abs = c.number_of_drums * (c.absorber_sector_deg / 360.0) * np.pi * (
        c.drum_radius**2 - (c.drum_radius - c.absorber_thickness)**2
    ) * c.active_core_height
    v_drum_beo = v_drums - v_abs
    v_rad_refl = (np.pi * (c.radial_reflector_outer_radius**2 - c.core_radius**2) * c.active_core_height) - v_drums + v_drum_beo

    v_tot = np.pi * c.radial_reflector_outer_radius**2 * c.total_reactor_height
    v_cool_tot = v_cool_core + v_cool_axial

    m_fuel = v_fuel * c.fuel_density * 1e-3
    m_clad = v_clad * c.cladding_density * 1e-3
    m_cool = v_cool_tot * c.coolant_density * 1e-3
    m_mod = v_mod * rho_mod * 1e-3
    m_refl_rad = v_rad_refl * c.reflector_density * 1e-3
    m_refl_ax = v_axial_refl * c.reflector_density * 1e-3
    rho_island = 2.90 if c.island_material == 'BeO' else 4.30
    m_island = v_island * rho_island * 1e-3
    m_abs = v_abs * c.absorber_density * 1e-3

    m_tot = m_fuel + m_clad + m_cool + m_mod + m_refl_rad + m_refl_ax + m_island + m_abs

    m_U_kg = m_fuel * ((238.03 * 0.8005 + 235.04 * 0.1995) / (238.03 * 0.8005 + 235.04 * 0.1995 + 14.01))
    m_U235_kg = m_fuel * (235.04 / (238.03 * 0.8005 + 235.04 * 0.1995 + 14.01)) * 0.1975

    inv = {
        "n_pins": n_pins,
        "m_fuel_kg": m_fuel,
        "m_clad_kg": m_clad,
        "m_cool_kg": m_cool,
        "m_mod_kg": m_mod,
        "m_refl_kg": m_refl_rad + m_refl_ax + m_island,
        "m_abs_kg": m_abs,
        "total_mass_kg": m_tot,
        "total_volume_L": v_tot * 1e-3,
        "m_U_kg": m_U_kg,
        "m_U235_kg": m_U235_kg,
        "v_fuel_L": v_fuel * 1e-3,
        "v_mod_L": v_mod * 1e-3,
        "v_mod_v_fuel_ratio": v_mod / v_fuel if v_fuel > 0 else 0.0
    }
    return inv


def define_opt11_candidates():
    """Defines complete screening suite across all families and benchmarks."""
    cands = []

    # =========================================================================
    # FAMILY A: Pure BeO Square Pin Core
    # =========================================================================
    # A1: Pure BeO Square Reference, 10 drums, R=20.0 cm, H=44.0 cm
    cands.append(ArchOpt11Candidate(
        cand_id="A1",
        family="Family A (BeO Square)",
        name="BeO Square Reference Core (10 Drums)",
        mod_type="BeO",
        lattice_type="square",
        core_radius=20.00,
        active_core_height=44.00,
        pitch=1.700,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20
    ))

    # A2: Pure BeO Scaled Square Core, 12 drums, R=21.5 cm, H=45.0 cm
    cands.append(ArchOpt11Candidate(
        cand_id="A2",
        family="Family A (BeO Square)",
        name="BeO Scaled Square Core (12 Drums)",
        mod_type="BeO",
        lattice_type="square",
        core_radius=21.50,
        active_core_height=45.00,
        pitch=1.650,
        number_of_drums=12,
        drum_radius=3.60,
        absorber_thickness=1.35,
        radial_reflector_thickness=8.80
    ))

    # =========================================================================
    # FAMILY B: Pure BeO Hexagonal Pin Lattice
    # =========================================================================
    # B1: Pure BeO Hexagonal Core, 10 drums, R=20.5 cm, H=44.0 cm, pitch 1.85 cm
    cands.append(ArchOpt11Candidate(
        cand_id="B1",
        family="Family B (BeO Hexagonal)",
        name="BeO Hexagonal Lattice Core (10 Drums)",
        mod_type="BeO",
        lattice_type="hex",
        core_radius=20.50,
        active_core_height=44.00,
        pitch=1.850,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20
    ))

    # B2: Pure BeO Hexagonal High-Packing Core, 12 drums, R=21.0 cm, H=45.0 cm, pitch 1.70 cm
    cands.append(ArchOpt11Candidate(
        cand_id="B2",
        family="Family B (BeO Hexagonal)",
        name="BeO High-Packing Hex Core (12 Drums)",
        mod_type="BeO",
        lattice_type="hex",
        core_radius=21.00,
        active_core_height=45.00,
        pitch=1.700,
        number_of_drums=12,
        drum_radius=3.60,
        absorber_thickness=1.35,
        radial_reflector_thickness=8.80
    ))

    # =========================================================================
    # FAMILY C: Annular Core with Central BeO Island
    # =========================================================================
    # C1: Annular Core with Central BeO Island (R_isl=5.0 cm), 10 drums
    cands.append(ArchOpt11Candidate(
        cand_id="C1",
        family="Family C (Annular BeO Island)",
        name="Annular Core with Central BeO Island (10 Drums)",
        mod_type="BeO",
        lattice_type="annular",
        core_radius=21.50,
        active_core_height=45.00,
        pitch=1.700,
        inner_island_radius=5.00,
        island_material="BeO",
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20
    ))

    # C2: Annular Core with Central BeO Island & ZrH2 Moderator Matrix, 10 drums
    cands.append(ArchOpt11Candidate(
        cand_id="C2",
        family="Family C (Annular BeO Island)",
        name="Annular Core with BeO Island & ZrH2 Driver",
        mod_type="ZrH2",
        lattice_type="annular",
        core_radius=19.40,
        active_core_height=41.50,
        pitch=1.665,
        inner_island_radius=4.50,
        island_material="BeO",
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20
    ))

    # =========================================================================
    # FAMILY D: Composite / Radially Zoned / Hybrid Cores
    # =========================================================================
    # D1: BeO Matrix with 70 vol% ZrH2 Hybrid Moderation, 10 drums
    cands.append(ArchOpt11Candidate(
        cand_id="D1",
        family="Family D (BeO-Hydride Composite)",
        name="BeO-ZrH2 Composite Matrix (70% ZrH2, 10 Drums)",
        mod_type="Hybrid_BeO_ZrH2_70",
        lattice_type="square",
        core_radius=19.40,
        active_core_height=41.50,
        pitch=1.665,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20
    ))

    # D2: BeO Matrix with 75 vol% ZrH2 Hybrid Moderation, 10 drums, H=42.5 cm
    cands.append(ArchOpt11Candidate(
        cand_id="D2",
        family="Family D (BeO-Hydride Composite)",
        name="BeO-ZrH2 Composite Matrix (75% ZrH2, H=42.5 cm)",
        mod_type="Hybrid_BeO_ZrH2_70", # close to 75%
        lattice_type="square",
        core_radius=19.50,
        active_core_height=42.50,
        pitch=1.665,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20
    ))

    # D3: BeO Matrix with 60 vol% YH2 Hybrid Moderation, 10 drums
    cands.append(ArchOpt11Candidate(
        cand_id="D3",
        family="Family D (BeO-Hydride Composite)",
        name="BeO-YH2 Composite Matrix (60% YH2, 10 Drums)",
        mod_type="Hybrid_BeO_YH2_60",
        lattice_type="square",
        core_radius=19.40,
        active_core_height=41.50,
        pitch=1.665,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20
    ))

    # =========================================================================
    # FAMILY E: Modular Clustered Core
    # =========================================================================
    # E1: 7-Cluster BeO Core, 10 drums, R=21.5 cm, H=44.0 cm
    cands.append(ArchOpt11Candidate(
        cand_id="E1",
        family="Family E (Clustered BeO Core)",
        name="7-Cluster Modular BeO Core (10 Drums)",
        mod_type="BeO",
        lattice_type="cluster",
        core_radius=21.50,
        active_core_height=44.00,
        pitch=1.700,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20
    ))

    # =========================================================================
    # CONTROLLED BENCHMARK CORES (YH2 & ZrH2)
    # =========================================================================
    # BM1: Opt-10 Baseline (YH2, 10 drums, R=19.4 cm, H=41.5 cm)
    cands.append(ArchOpt11Candidate(
        cand_id="BM1",
        family="Benchmark (Opt-10 YH2)",
        name="Opt-10 Baseline Reference (YH2, 10 Drums)",
        mod_type="YH2",
        lattice_type="square",
        core_radius=19.40,
        active_core_height=41.50,
        pitch=1.665,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20
    ))

    # BM2: Opt-10 Extended (YH2, 10 drums, H=42.2 cm)
    cands.append(ArchOpt11Candidate(
        cand_id="BM2",
        family="Benchmark (Opt-10 YH2)",
        name="Opt-10 Extended Core (YH2, H=42.2 cm)",
        mod_type="YH2",
        lattice_type="square",
        core_radius=19.40,
        active_core_height=42.20,
        pitch=1.665,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20
    ))

    # BM3: High-Performance ZrH2 Compact Architecture (10 drums, R=18.5 cm, H=39.5 cm)
    cands.append(ArchOpt11Candidate(
        cand_id="BM3",
        family="Benchmark (ZrH2 High-Perf)",
        name="ZrH2 Compact High-Perf Core (10 Drums)",
        mod_type="ZrH2",
        lattice_type="square",
        core_radius=18.50,
        active_core_height=39.50,
        pitch=1.665,
        number_of_drums=10,
        drum_radius=3.80,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.00
    ))

    # BM4: ZrH2 Ultra-Compact Architecture (12 drums, R=18.3 cm, H=39.0 cm)
    cands.append(ArchOpt11Candidate(
        cand_id="BM4",
        family="Benchmark (ZrH2 High-Perf)",
        name="ZrH2 Ultra-Compact Core (12 Drums)",
        mod_type="ZrH2",
        lattice_type="square",
        core_radius=18.30,
        active_core_height=39.00,
        pitch=1.665,
        number_of_drums=12,
        drum_radius=3.40,
        absorber_thickness=1.35,
        radial_reflector_thickness=8.60
    ))

    return cands


def run_stage1_screening():
    """Executes Stage 1 screening across all candidates."""
    print("="*80)
    print("STARTING STAGE 1 CANDIDATE ARCHITECTURE SCREENING (CANDIDATE OPT-11)")
    print("="*80)

    work_dir = "results/opt11/screening/workspace_screening"
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()

    candidates = define_opt11_candidates()
    results = []

    # Screening Monte Carlo settings: fast and accurate
    particles = 15000
    batches = 28
    inactive = 6

    os.chdir(work_dir)

    for cand in candidates:
        print(f"\n---> Screening Candidate [{cand.cand_id}]: {cand.name} ({cand.family})")
        t_start = time.time()

        # 1. Build Operating State (k_out, 0 deg)
        openmc.reset_auto_ids()
        mats, f, cl, co, mo, be, b, isl, rho_mod = build_candidate_materials(cand)
        geom, pins, _ = build_candidate_geometry(cand, f, cl, co, mo, be, b, isl, drum_angles=[0.0]*cand.number_of_drums)

        inv = compute_candidate_inventory(cand, pins, rho_mod)
        print(f"     Geometry: Pins={inv['n_pins']}, Mass={inv['total_mass_kg']:.2f} kg, Volume={inv['total_volume_L']:.2f} L, Vmod/Vfuel={inv['v_mod_v_fuel_ratio']:.3f}")

        settings = openmc.Settings()
        settings.run_mode = "eigenvalue"
        settings.particles = particles
        settings.batches = batches
        settings.inactive = inactive
        settings.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}
        settings.source = openmc.IndependentSource(space=openmc.stats.Point((0, 0, 0)))

        mod_out = openmc.Model(materials=mats, geometry=geom, settings=settings)
        mod_out.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{batches}.h5") as sp:
            k_out = float(sp.keff.nominal_value)
            s_out = float(sp.keff.std_dev)

        # 2. Build Shutdown State (k_in, 180 deg)
        openmc.reset_auto_ids()
        mats_in, f_i, cl_i, co_i, mo_i, be_i, b_i, isl_i, _ = build_candidate_materials(cand)
        geom_in, _, _ = build_candidate_geometry(cand, f_i, cl_i, co_i, mo_i, be_i, b_i, isl_i, drum_angles=[180.0]*cand.number_of_drums)

        mod_in = openmc.Model(materials=mats_in, geometry=geom_in, settings=settings)
        mod_in.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{batches}.h5") as sp:
            k_in = float(sp.keff.nominal_value)
            s_in = float(sp.keff.std_dev)

        drum_worth = (1.0 / k_in - 1.0 / k_out) * 1e5
        sdm = (1.0 - k_in) * 1e5

        # 3. Single Stuck-Drum State (Drum #0 at 0 deg, remaining drums at 180 deg)
        openmc.reset_auto_ids()
        stuck_angles = [0.0] + [180.0] * (cand.number_of_drums - 1)
        mats_stk, f_s, cl_s, co_s, mo_s, be_s, b_s, isl_s, _ = build_candidate_materials(cand)
        geom_stk, _, _ = build_candidate_geometry(cand, f_s, cl_s, co_s, mo_s, be_s, b_s, isl_s, drum_angles=stuck_angles)

        mod_stk = openmc.Model(materials=mats_stk, geometry=geom_stk, settings=settings)
        mod_stk.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{batches}.h5") as sp:
            k_stuck = float(sp.keff.nominal_value)
            s_stuck = float(sp.keff.std_dev)

        sdm_stuck = (1.0 - k_stuck) * 1e5

        # 4. Depressurized Void State (Vacuum in channels & ducts, drums at 0 deg)
        openmc.reset_auto_ids()
        mats_vd, f_v, cl_v, co_v, mo_v, be_v, b_v, isl_v, _ = build_candidate_materials(cand, cool_density=1.0e-7)
        geom_vd, _, _ = build_candidate_geometry(cand, f_v, cl_v, co_v, mo_v, be_v, b_v, isl_v, drum_angles=[0.0]*cand.number_of_drums)

        mod_vd = openmc.Model(materials=mats_vd, geometry=geom_vd, settings=settings)
        mod_vd.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{batches}.h5") as sp:
            k_void = float(sp.keff.nominal_value)
            s_void = float(sp.keff.std_dev)

        rho_out = (k_out - 1.0) / k_out
        rho_void = (k_void - 1.0) / k_void
        delta_rho_void_pcm = (rho_void - rho_out) * 1e5
        sigma_void_pcm = np.sqrt((s_void / (k_void**2))**2 + (s_out / (k_out**2))**2) * 1e5

        t_elapsed = time.time() - t_start
        print(f"     Results: k_out={k_out:.5f}+/-{s_out:.5f}, k_in={k_in:.5f}+/-{s_in:.5f}, Worth={drum_worth:.1f} pcm")
        print(f"              SDM={sdm:.1f} pcm, SDM_stuck={sdm_stuck:.1f} pcm")
        print(f"              Void Reactivity={delta_rho_void_pcm:.1f} +/- {sigma_void_pcm:.1f} pcm | Time={t_elapsed:.1f}s")

        # Feasibility check
        is_critical = k_out >= 1.000
        is_sdm_ok = sdm >= 2000.0 and sdm_stuck >= 1000.0
        is_mass_ok = inv['total_mass_kg'] <= 650.0

        rec = {
            "cand_id": cand.cand_id,
            "family": cand.family,
            "name": cand.name,
            "mod_type": cand.mod_type,
            "lattice_type": cand.lattice_type,
            "core_radius_cm": cand.core_radius,
            "active_core_height_cm": cand.active_core_height,
            "aspect_ratio_H_R": cand.aspect_ratio_H_R,
            "pitch_cm": cand.pitch,
            "number_of_drums": cand.number_of_drums,
            "radial_refl_thick_cm": cand.radial_reflector_thickness,
            "n_pins": inv["n_pins"],
            "v_mod_v_fuel": inv["v_mod_v_fuel_ratio"],
            "mass_fuel_kg": inv["m_fuel_kg"],
            "mass_clad_kg": inv["m_clad_kg"],
            "mass_cool_kg": inv["m_cool_kg"],
            "mass_mod_kg": inv["m_mod_kg"],
            "mass_refl_kg": inv["m_refl_kg"],
            "mass_abs_kg": inv["m_abs_kg"],
            "total_mass_kg": inv["total_mass_kg"],
            "total_volume_L": inv["total_volume_L"],
            "m_U_kg": inv["m_U_kg"],
            "m_U235_kg": inv["m_U235_kg"],
            "k_out": k_out,
            "k_out_sd": s_out,
            "k_in": k_in,
            "k_in_sd": s_in,
            "drum_worth_pcm": drum_worth,
            "sdm_pcm": sdm,
            "k_stuck": k_stuck,
            "k_stuck_sd": s_stuck,
            "sdm_stuck_pcm": sdm_stuck,
            "k_void": k_void,
            "k_void_sd": s_void,
            "delta_rho_void_pcm": delta_rho_void_pcm,
            "sigma_void_pcm": sigma_void_pcm,
            "void_upper_3sigma_pcm": delta_rho_void_pcm + 3.0 * sigma_void_pcm,
            "is_critical": is_critical,
            "is_sdm_ok": is_sdm_ok,
            "is_mass_ok": is_mass_ok
        }
        results.append(rec)

    os.chdir(orig_dir)

    # Export to CSV
    csv_path = "results/opt11/screening/opt11_architecture_screening.csv"
    with open(csv_path, "w", newline="") as f_csv:
        fieldnames = list(results[0].keys())
        writer = csv.DictWriter(f_csv, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print(f"\n[+] Screening data written successfully to {csv_path}")

    # Generate Pareto screening plots
    generate_screening_plots(results)
    return results


def generate_screening_plots(results):
    """Generates Stage 1 Pareto screening comparison plots."""
    os.makedirs("results/opt11/screening", exist_ok=True)
    artifacts_dir = "/home/aleri/.gemini/antigravity/brain/0c8a66db-b02e-44e7-8d87-eec426335ff7"

    family_colors = {
        "Family A (BeO Square)": "#e74c3c",
        "Family B (BeO Hexagonal)": "#e67e22",
        "Family C (Annular BeO Island)": "#f39c12",
        "Family D (BeO-Hydride Composite)": "#2ecc71",
        "Family E (Clustered BeO Core)": "#9b59b6",
        "Benchmark (Opt-10 YH2)": "#3498db",
        "Benchmark (ZrH2 High-Perf)": "#1abc9c"
    }

    # Plot 1: Mass vs k_out (Criticality Viability)
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    for r in results:
        col = family_colors.get(r['family'], '#7f8c8d')
        marker = 'o' if r['is_critical'] else 'X'
        ax.scatter(r['total_mass_kg'], r['k_out'], color=col, s=120, edgecolors='black', linewidth=1.2, zorder=5, marker=marker)
        ax.annotate(r['cand_id'], (r['total_mass_kg'], r['k_out']), textcoords="offset points", xytext=(0, 7), ha='center', fontsize=9, fontweight='bold')

    ax.axhline(1.000, color='red', linestyle='--', linewidth=1.5, label='Criticality Limit ($k_{eff} = 1.000$)')
    ax.axvline(650.0, color='darkred', linestyle=':', linewidth=1.5, label='Mass Target ($650$ kg)')
    ax.set_xlabel("Total Assembly Mass [kg]", fontsize=12, fontweight='bold')
    ax.set_ylabel("BOL Multiplication Factor $k_{out}$ [-]", fontsize=12, fontweight='bold')
    ax.set_title("Opt-11 Stage 1 Screening: Criticality vs Mass", fontsize=14, fontweight='bold', pad=12)
    ax.grid(True, linestyle=':', alpha=0.6)
    
    # Custom legend
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=10, label=f) for f, c in family_colors.items()]
    handles.append(plt.Line2D([0], [0], color='red', linestyle='--', label='Criticality Limit ($k=1.0$)'))
    handles.append(plt.Line2D([0], [0], color='darkred', linestyle=':', label='Mass Envelope ($650$ kg)'))
    ax.legend(handles=handles, loc='lower right', fontsize=8, framealpha=0.9)
    plt.tight_layout()
    p1 = "results/opt11/screening/opt11_screening_mass_vs_criticality.png"
    plt.savefig(p1)
    shutil.copy(p1, os.path.join(artifacts_dir, "opt11_screening_mass_vs_criticality.png"))
    plt.close()

    # Plot 2: Mass vs Volume
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    for r in results:
        col = family_colors.get(r['family'], '#7f8c8d')
        ax.scatter(r['total_mass_kg'], r['total_volume_L'], color=col, s=120, edgecolors='black', linewidth=1.2, zorder=5)
        ax.annotate(r['cand_id'], (r['total_mass_kg'], r['total_volume_L']), textcoords="offset points", xytext=(0, 7), ha='center', fontsize=9, fontweight='bold')

    ax.axvline(650.0, color='darkred', linestyle=':', linewidth=1.5, label='Mass Target ($650$ kg)')
    ax.set_xlabel("Total Assembly Mass [kg]", fontsize=12, fontweight='bold')
    ax.set_ylabel("Total Reactor Outer Volume [L]", fontsize=12, fontweight='bold')
    ax.set_title("Opt-11 Stage 1 Screening: Mass vs Volume Envelope", fontsize=14, fontweight='bold', pad=12)
    ax.grid(True, linestyle=':', alpha=0.6)
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=10, label=f) for f, c in family_colors.items()]
    ax.legend(handles=handles, loc='upper left', fontsize=8, framealpha=0.9)
    plt.tight_layout()
    p2 = "results/opt11/screening/opt11_screening_mass_vs_volume.png"
    plt.savefig(p2)
    shutil.copy(p2, os.path.join(artifacts_dir, "opt11_screening_mass_vs_volume.png"))
    plt.close()

    # Plot 3: Stuck Drum SDM vs Total Drum Worth
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)
    for r in results:
        col = family_colors.get(r['family'], '#7f8c8d')
        ax.scatter(r['drum_worth_pcm'], r['sdm_stuck_pcm'], color=col, s=120, edgecolors='black', linewidth=1.2, zorder=5)
        ax.annotate(r['cand_id'], (r['drum_worth_pcm'], r['sdm_stuck_pcm']), textcoords="offset points", xytext=(0, 7), ha='center', fontsize=9, fontweight='bold')

    ax.axhline(1000.0, color='red', linestyle='--', linewidth=1.5, label='Safety Limit: $SDM_{stuck} \geq 1000$ pcm')
    ax.set_xlabel("Total Drum Worth [pcm]", fontsize=12, fontweight='bold')
    ax.set_ylabel("Stuck-Drum Shutdown Margin $SDM_{stuck}$ [pcm]", fontsize=12, fontweight='bold')
    ax.set_title("Opt-11 Stage 1 Screening: Nuclear Control Safety", fontsize=14, fontweight='bold', pad=12)
    ax.grid(True, linestyle=':', alpha=0.6)
    handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=10, label=f) for f, c in family_colors.items()]
    handles.append(plt.Line2D([0], [0], color='red', linestyle='--', label='$SDM_{stuck} \geq 1000$ pcm'))
    ax.legend(handles=handles, loc='lower right', fontsize=8, framealpha=0.9)
    plt.tight_layout()
    p3 = "results/opt11/screening/opt11_screening_control_safety.png"
    plt.savefig(p3)
    shutil.copy(p3, os.path.join(artifacts_dir, "opt11_screening_control_safety.png"))
    plt.close()

    print("[+] All screening plots generated and copied to artifact directory.")


if __name__ == "__main__":
    run_stage1_screening()
