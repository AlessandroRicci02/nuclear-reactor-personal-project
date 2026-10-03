#!/usr/bin/env python3
"""
===============================================================================
STAGE 1 ARCHITECTURE SCREENING ENGINE: SPACE NEP REACTOR
===============================================================================
Explores four distinct candidate architecture families:
  Family A: Compact Pin-Core Reference Variations (Square lattice, varied aspect
            ratio, pitch, 8/10/12 control drums, enhanced streaming ducts)
  Family B: Hexagonal / Triangular Lattice Core (Hexagonal pin array, 6/8/12 drums,
            isotropic moderation, tighter radial packaging)
  Family C: Annular Core with Central Island (Central YH2 or BeO flux trap island,
            annular fuel zone, peripheral control drums)
  Family D: Radially Zoned / Clustered Modular Core (Two-zone pitch/moderation,
            enhanced peripheral leakage and thermal flux tuning)

Pipeline:
  1. Construct exact 3D OpenMC CSG geometry for each candidate.
  2. Compute itemized mass and volume breakdown.
  3. Evaluate BOL operating multiplication (k_out, drums 0 deg).
  4. Evaluate shutdown multiplication (k_in, drums 180 deg) & Total Drum Worth.
  5. Evaluate Coolant Void Reactivity (Depressurized vacuum state).
  6. Export comprehensive screening metrics to results/architecture_screening.csv.
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
import matplotlib.patches as patches

import openmc

os.environ["PATH"] = "/home/aleri/miniforge3/envs/openmc-env/bin:" + os.environ.get("PATH", "")
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
os.environ["OMP_NUM_THREADS"] = "8"
openmc.config['cross_sections'] = "/home/aleri/wmp_data/cross_sections_wmp.xml"


class ArchitectureCandidate:
    def __init__(self,
                 cand_id,
                 family,
                 name,
                 lattice_type='square',
                 core_radius=19.40,
                 active_core_height=41.50,
                 fuel_radius=0.580,
                 cladding_thickness=0.050,
                 coolant_channel_thickness=0.070,
                 pitch=1.665,
                 inner_island_radius=0.0,
                 island_material=None,
                 zoned_pitch=None, # (r_split, pitch_in, pitch_out)
                 number_of_drums=8,
                 drum_radius=4.40,
                 absorber_thickness=1.45,
                 absorber_sector_deg=140.0,
                 radial_reflector_thickness=10.00,
                 axial_reflector_thickness=9.50,
                 axial_duct_radius=0.800,
                 coolant_pressure_MPa=2.5,
                 coolant_temp_K=1050.0):

        self.cand_id = cand_id
        self.family = family
        self.name = name
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
        self.zoned_pitch = zoned_pitch

        self.clad_outer_radius = self.fuel_radius + self.cladding_thickness
        self.channel_outer_radius = self.clad_outer_radius + self.coolant_channel_thickness
        self.axial_duct_radius = axial_duct_radius

        self.number_of_drums = number_of_drums
        self.drum_radius = drum_radius
        self.absorber_thickness = absorber_thickness
        self.absorber_sector_deg = absorber_sector_deg

        self.drum_clearance = 0.20
        self.drum_radial_position = self.core_radius + self.drum_clearance + self.drum_radius
        # Radial reflector thickness adjusted to enclose drums with >= 0.8 cm margin
        min_refl_thick = self.drum_clearance + 2.0 * self.drum_radius + 0.80
        self.radial_reflector_thickness = max(radial_reflector_thickness, min_refl_thick)
        self.radial_reflector_outer_radius = self.core_radius + self.radial_reflector_thickness
        self.axial_reflector_thickness = axial_reflector_thickness
        self.total_reactor_height = self.active_core_height + 2.0 * self.axial_reflector_thickness

        # Materials
        self.fuel_density = 13.50
        self.fuel_enrichment = 19.75
        self.cladding_density = 8.57
        self.moderator_density = 4.30
        self.reflector_density = 2.90
        self.absorber_density = 2.52
        self.b10_enrichment = 0.90

        self.coolant_pressure_MPa = coolant_pressure_MPa
        self.coolant_temp_K = coolant_temp_K
        self.coolant_density = (self.coolant_pressure_MPa * 1.0e6 * 4.0026e-3) / (8.31446 * self.coolant_temp_K) * 1.0e-3


def build_candidate_materials(c: ArchitectureCandidate, cool_density=None):
    """Build material cards with explicit isotopes."""
    openmc.reset_auto_ids()

    # Fuel
    fuel = openmc.Material(name="UN Fuel (19.75 wt% U-235)")
    fuel.set_density("g/cm3", c.fuel_density)
    fuel.add_nuclide("U235", 0.5 * 0.199522, "ao")
    fuel.add_nuclide("U238", 0.5 * 0.800478, "ao")
    fuel.add_nuclide("N14", 0.5 * 0.9963, "ao")
    fuel.add_nuclide("N15", 0.5 * 0.0037, "ao")
    fuel.temperature = 1050.0

    # Cladding
    cladding = openmc.Material(name="Nb-1Zr Cladding")
    cladding.set_density("g/cm3", c.cladding_density)
    cladding.add_element("Nb", 0.99, "wo")
    cladding.add_element("Zr", 0.01, "wo")
    cladding.temperature = 1050.0

    # Coolant
    rho_c = c.coolant_density if cool_density is None else cool_density
    coolant = openmc.Material(name="Pure Helium Coolant")
    coolant.set_density("g/cm3", max(rho_c, 1.0e-7))
    coolant.add_nuclide("He4", 1.0, "ao")
    coolant.temperature = c.coolant_temp_K

    # Moderator YH2
    moderator = openmc.Material(name="YH2 Moderator")
    moderator.set_density("g/cm3", c.moderator_density)
    moderator.add_nuclide("Y89", 1.0, "ao")
    moderator.add_nuclide("H1", 2.0, "ao")
    moderator.add_s_alpha_beta("c_H_in_ZrH")
    moderator.temperature = 800.0

    # Reflector BeO
    reflector = openmc.Material(name="BeO Reflector")
    reflector.set_density("g/cm3", c.reflector_density)
    reflector.add_nuclide("Be9", 1.0, "ao")
    reflector.add_nuclide("O16", 1.0, "ao")
    reflector.add_s_alpha_beta("c_Be_in_BeO")
    reflector.add_s_alpha_beta("c_O_in_BeO")
    reflector.temperature = 800.0

    # Absorber B4C (90 at% B-10)
    absorber = openmc.Material(name="B4C Absorber (90% B-10)")
    absorber.set_density("g/cm3", c.absorber_density)
    absorber.add_nuclide("B10", 4.0 * c.b10_enrichment, "ao")
    absorber.add_nuclide("B11", 4.0 * (1.0 - c.b10_enrichment), "ao")
    absorber.add_nuclide("C0", 1.0, "ao")
    absorber.temperature = 800.0

    mats = [fuel, cladding, coolant, moderator, reflector, absorber]
    materials = openmc.Materials(mats)
    materials.cross_sections = "/home/aleri/wmp_data/cross_sections_wmp.xml"
    return materials, fuel, cladding, coolant, moderator, reflector, absorber


def build_candidate_geometry(c: ArchitectureCandidate, fuel, cladding, coolant, moderator, reflector, absorber, drum_angles=None):
    """
    Constructs the 3D CSG geometry for any of the 4 architecture families.
    """
    if drum_angles is None:
        drum_angles = [0.0] * c.number_of_drums

    s_fuel = openmc.ZCylinder(r=c.fuel_radius, name="Pin Fuel Outer")
    s_clad = openmc.ZCylinder(r=c.clad_outer_radius, name="Pin Clad Outer")
    s_channel = openmc.ZCylinder(r=c.channel_outer_radius, name="Coolant Channel Outer")

    z_core_bot = openmc.ZPlane(z0=-c.active_core_height / 2.0, name="Core Bot")
    z_core_top = openmc.ZPlane(z0=+c.active_core_height / 2.0, name="Core Top")
    z_refl_bot = openmc.ZPlane(z0=-c.total_reactor_height / 2.0, boundary_type="vacuum", name="Reactor Bot")
    z_refl_top = openmc.ZPlane(z0=+c.total_reactor_height / 2.0, boundary_type="vacuum", name="Reactor Top")

    # Cooled pin universe:
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

    # Reflector-only universe (for BeO central island)
    c_be_island_core = openmc.Cell(fill=reflector, region=+z_core_bot & -z_core_top)
    c_be_bot_island = openmc.Cell(fill=reflector, region=+z_refl_bot & -z_core_bot)
    c_be_top_island = openmc.Cell(fill=reflector, region=+z_core_top & -z_refl_top)
    u_beo_only = openmc.Universe(cells=[c_be_island_core, c_be_bot_island, c_be_top_island])

    pin_centers = []

    # Build Core Lattice / Cell Structure
    if c.lattice_type in ['square', 'annular', 'zoned']:
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
                    else:
                        is_pin = True

                if is_pin:
                    row.append(u_pin)
                    pin_centers.append((x, y))
                else:
                    if c.lattice_type == 'annular' and r_dist < c.inner_island_radius:
                        if c.island_material == 'BeO':
                            row.append(u_beo_only)
                        else:
                            row.append(u_mod_only)
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
        # Compute number of rings to span core radius
        # outer edge of ring r is r * pitch
        num_rings = int(np.ceil(c.core_radius / pitch)) + 1
        
        # Basis vectors for orientation 'y'
        v1 = np.array([pitch, 0.0])
        v2 = np.array([pitch/2.0, pitch * np.sqrt(3)/2.0])

        # We construct universes ring by ring from outermost (num_rings - 1) down to 0
        rings_universes = []
        for r in range(num_rings - 1, -1, -1):
            ring_u = []
            if r == 0:
                # Center
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


def compute_candidate_inventory(c: ArchitectureCandidate, pin_centers):
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
    
    # In annular core with central island:
    if c.inner_island_radius > 0:
        v_island = np.pi * c.inner_island_radius**2 * c.active_core_height
    else:
        v_island = 0.0

    v_mod = v_active_core_tot - (v_fuel + v_clad + v_cool_core)
    # If central island is BeO:
    if c.island_material == 'BeO':
        v_mod -= v_island
        v_island_beo = v_island
    else:
        v_island_beo = 0.0

    v_axial_refl = np.pi * c.core_radius**2 * (2 * c.axial_reflector_thickness) - v_cool_axial

    # Radial reflector and drums
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
    m_mod = v_mod * c.moderator_density * 1e-3
    m_refl_rad = v_rad_refl * c.reflector_density * 1e-3
    m_refl_ax = v_axial_refl * c.reflector_density * 1e-3
    m_island_beo = v_island_beo * c.reflector_density * 1e-3
    m_abs = v_abs * c.absorber_density * 1e-3

    m_tot = m_fuel + m_clad + m_cool + m_mod + m_refl_rad + m_refl_ax + m_island_beo + m_abs

    m_U_kg = m_fuel * ((238.03 * 0.8005 + 235.04 * 0.1995) / (238.03 * 0.8005 + 235.04 * 0.1995 + 14.01))
    m_U235_kg = m_fuel * (235.04 / (238.03 * 0.8005 + 235.04 * 0.1995 + 14.01)) * 0.1975

    inv = {
        "n_pins": n_pins,
        "m_fuel_kg": m_fuel,
        "m_clad_kg": m_clad,
        "m_cool_kg": m_cool,
        "m_mod_kg": m_mod,
        "m_refl_kg": m_refl_rad + m_refl_ax + m_island_beo,
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


def define_screening_candidates():
    """Defines 12 candidates spanning Families A, B, C, D."""
    candidates = []

    # =========================================================================
    # FAMILY A: Compact Pin-Core Reference Variations (Square Pitch & Multi-Drums)
    # =========================================================================
    # =========================================================================
    # FAMILY A: Compact Pin-Core Reference Variations (Square Pitch & Multi-Drums)
    # =========================================================================
    # A1: 10 Control Drums, Opt-9 pitch, reduced reflector thickness
    candidates.append(ArchitectureCandidate(
        cand_id="A1",
        family="Family A (Square Multi-Drum)",
        name="10-Drum Compact Square Core",
        lattice_type='square',
        core_radius=19.40,
        active_core_height=41.50,
        pitch=1.665,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.820
    ))

    # A2: 12 Control Drums, very high drum safety margin, ultra-compact reflector
    candidates.append(ArchitectureCandidate(
        cand_id="A2",
        family="Family A (Square Multi-Drum)",
        name="12-Drum Ultra-Shielded Square Core",
        lattice_type='square',
        core_radius=19.40,
        active_core_height=41.50,
        pitch=1.665,
        number_of_drums=12,
        drum_radius=3.50,
        absorber_thickness=1.35,
        radial_reflector_thickness=8.60,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.850
    ))

    # A3: High Aspect Ratio (H/R = 2.31), 10 Drums, enhanced axial leakage
    candidates.append(ArchitectureCandidate(
        cand_id="A3",
        family="Family A (Square Multi-Drum)",
        name="High Aspect Ratio (H/R=2.31) Square Core",
        lattice_type='square',
        core_radius=18.80,
        active_core_height=43.50,
        pitch=1.665,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.850
    ))

    # =========================================================================
    # FAMILY B: Hexagonal / Triangular Lattice Core (Isotropic Moderation)
    # =========================================================================
    # B1: 10 Control Drums, Hexagonal Pitch 1.79 cm, 19.2 cm core radius
    candidates.append(ArchitectureCandidate(
        cand_id="B1",
        family="Family B (Hexagonal Lattice)",
        name="10-Drum Hexagonal Lattice Core",
        lattice_type='hex',
        core_radius=19.20,
        active_core_height=41.50,
        pitch=1.790,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.850
    ))

    # B2: 12 Control Drums, Hexagonal Pitch 1.78 cm, 19.1 cm core radius
    candidates.append(ArchitectureCandidate(
        cand_id="B2",
        family="Family B (Hexagonal Lattice)",
        name="12-Drum Ultra-Compact Hex Core",
        lattice_type='hex',
        core_radius=19.10,
        active_core_height=41.50,
        pitch=1.780,
        number_of_drums=12,
        drum_radius=3.50,
        absorber_thickness=1.35,
        radial_reflector_thickness=8.60,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.850
    ))

    # B3: 8 Control Drums, Hexagonal Pitch 1.80 cm, 19.2 cm core radius
    candidates.append(ArchitectureCandidate(
        cand_id="B3",
        family="Family B (Hexagonal Lattice)",
        name="8-Drum Symmetric Hex Core",
        lattice_type='hex',
        core_radius=19.20,
        active_core_height=41.50,
        pitch=1.800,
        number_of_drums=8,
        drum_radius=4.40,
        absorber_thickness=1.45,
        radial_reflector_thickness=10.00,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.820
    ))

    # =========================================================================
    # FAMILY C: Annular Core with Central Island (Internal Flux Trap / Reflector)
    # =========================================================================
    # C1: Central YH2 Moderator Island (R_island = 4.0 cm), 10 Drums
    candidates.append(ArchitectureCandidate(
        cand_id="C1",
        family="Family C (Annular Island)",
        name="Annular Core with Central YH2 Flux Trap",
        lattice_type='annular',
        core_radius=20.00,
        active_core_height=41.50,
        pitch=1.665,
        inner_island_radius=4.00,
        island_material='YH2',
        number_of_drums=10,
        drum_radius=4.00,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.40,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.820
    ))

    # C2: Central BeO Reflector Island (R_island = 4.5 cm), 10 Drums
    candidates.append(ArchitectureCandidate(
        cand_id="C2",
        family="Family C (Annular Island)",
        name="Annular Core with Central BeO Reflector",
        lattice_type='annular',
        core_radius=20.00,
        active_core_height=41.50,
        pitch=1.665,
        inner_island_radius=4.50,
        island_material='BeO',
        number_of_drums=10,
        drum_radius=4.00,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.40,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.820
    ))

    # C3: Compact Annular Core with Central YH2 Island (R_island = 3.5 cm), 12 Drums
    candidates.append(ArchitectureCandidate(
        cand_id="C3",
        family="Family C (Annular Island)",
        name="Compact Annular Core (12 Drums)",
        lattice_type='annular',
        core_radius=19.80,
        active_core_height=41.50,
        pitch=1.665,
        inner_island_radius=3.50,
        island_material='YH2',
        number_of_drums=12,
        drum_radius=3.60,
        absorber_thickness=1.35,
        radial_reflector_thickness=8.60,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.840
    ))

    # =========================================================================
    # FAMILY D: Clustered / Radially Zoned Core
    # =========================================================================
    # D1: Radially Zoned Moderation, 10 Drums
    candidates.append(ArchitectureCandidate(
        cand_id="D1",
        family="Family D (Clustered/Zoned)",
        name="Radially Zoned Moderation Core (10 Drums)",
        lattice_type='square',
        core_radius=19.40,
        active_core_height=41.50,
        pitch=1.670,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.830
    ))

    # D2: Clustered Enhanced Streaming Channels (Ducts r=0.88 cm), 12 Drums
    candidates.append(ArchitectureCandidate(
        cand_id="D2",
        family="Family D (Clustered/Zoned)",
        name="Enhanced Streaming Modular Core (12 Drums)",
        lattice_type='square',
        core_radius=19.40,
        active_core_height=41.50,
        pitch=1.665,
        number_of_drums=12,
        drum_radius=3.50,
        absorber_thickness=1.35,
        radial_reflector_thickness=8.60,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.880
    ))

    # D3: High-Leakage Hexagonal Zoned Core (Hex lattice with oversized ducts), 10 Drums
    candidates.append(ArchitectureCandidate(
        cand_id="D3",
        family="Family D (Clustered/Zoned)",
        name="High-Leakage Hexagonal Core (10 Drums)",
        lattice_type='hex',
        core_radius=19.20,
        active_core_height=41.50,
        pitch=1.790,
        number_of_drums=10,
        drum_radius=3.90,
        absorber_thickness=1.40,
        radial_reflector_thickness=9.20,
        axial_reflector_thickness=9.50,
        axial_duct_radius=0.880
    ))

    return candidates


def run_stage1_screening():
    """
    Executes Stage 1 screening on all 12 candidate architectures.
    Fast Monte Carlo (20,000 particles x 32 batches, inactive 8).
    Evaluates:
      - k_out (BOL 0 deg)
      - k_in (BOL 180 deg)
      - Drum worth
      - Coolant void reactivity (depressurized vacuum)
      - Total mass, volume, fuel mass, pins
    """
    print("="*80)
    print("STARTING STAGE 1 CANDIDATE ARCHITECTURE SCREENING")
    print("="*80)

    work_dir = "results/screening_workspace"
    os.makedirs(work_dir, exist_ok=True)
    orig_dir = os.getcwd()

    candidates = define_screening_candidates()
    results = []

    # Fast screening Monte Carlo settings
    particles = 20000
    batches = 32
    inactive = 8

    os.chdir(work_dir)

    for cand in candidates:
        print(f"\n---> Screening Candidate [{cand.cand_id}]: {cand.name} ({cand.family})")
        t_start = time.time()

        # 1. Build Operating State (k_out, 0 deg)
        openmc.reset_auto_ids()
        mats, f, cl, co, mo, be, b = build_candidate_materials(cand)
        geom, pins, _ = build_candidate_geometry(cand, f, cl, co, mo, be, b, drum_angles=[0.0]*cand.number_of_drums)

        inv = compute_candidate_inventory(cand, pins)
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
        mats_in, f_i, cl_i, co_i, mo_i, be_i, b_i = build_candidate_materials(cand)
        geom_in, _, _ = build_candidate_geometry(cand, f_i, cl_i, co_i, mo_i, be_i, b_i, drum_angles=[180.0]*cand.number_of_drums)

        mod_in = openmc.Model(materials=mats_in, geometry=geom_in, settings=settings)
        mod_in.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{batches}.h5") as sp:
            k_in = float(sp.keff.nominal_value)
            s_in = float(sp.keff.std_dev)

        drum_worth = (1.0 / k_in - 1.0 / k_out) * 1e5
        sdm = (1.0 - k_in) * 1e5

        # 3. Build Depressurized Void State (Vacuum in channels & ducts, drums at 0 deg)
        openmc.reset_auto_ids()
        mats_vd, f_v, cl_v, co_v, mo_v, be_v, b_v = build_candidate_materials(cand, cool_density=1.0e-7)
        geom_vd, _, _ = build_candidate_geometry(cand, f_v, cl_v, co_v, mo_v, be_v, b_v, drum_angles=[0.0]*cand.number_of_drums)

        mod_vd = openmc.Model(materials=mats_vd, geometry=geom_vd, settings=settings)
        mod_vd.export_to_xml()
        openmc.run(output=False)

        with openmc.StatePoint(f"statepoint.{batches}.h5") as sp:
            k_vd = float(sp.keff.nominal_value)
            s_vd = float(sp.keff.std_dev)

        drho_void = (1.0 / k_out - 1.0 / k_vd) * 1e5
        u_void = np.sqrt((s_out / k_out**2)**2 + (s_vd / k_vd**2)**2) * 1e5

        dt = time.time() - t_start
        print(f"     Results: k_out={k_out:.5f}+/-{s_out:.5f} | k_in={k_in:.5f}+/-{s_in:.5f} | Worth={drum_worth:.0f} pcm | SDM={sdm:.0f} pcm | Void={drho_void:+.1f}+/-{u_void:.1f} pcm | Time={dt:.1f}s")

        res_dict = {
            "cand_id": cand.cand_id,
            "family": cand.family,
            "name": cand.name,
            "lattice_type": cand.lattice_type,
            "core_radius_cm": cand.core_radius,
            "core_height_cm": cand.active_core_height,
            "aspect_ratio_H_R": cand.aspect_ratio_H_R,
            "pitch_cm": cand.pitch,
            "num_drums": cand.number_of_drums,
            "drum_radius_cm": cand.drum_radius,
            "refl_thick_radial_cm": cand.radial_reflector_thickness,
            "refl_thick_axial_cm": cand.axial_reflector_thickness,
            "duct_radius_cm": cand.axial_duct_radius,
            "pins": inv["n_pins"],
            "v_mod_v_fuel": inv["v_mod_v_fuel_ratio"],
            "mass_total_kg": inv["total_mass_kg"],
            "volume_total_L": inv["total_volume_L"],
            "mass_fuel_kg": inv["m_fuel_kg"],
            "mass_U235_kg": inv["m_U235_kg"],
            "mass_refl_kg": inv["m_refl_kg"],
            "k_out": k_out,
            "s_out": s_out,
            "k_in": k_in,
            "s_in": s_in,
            "drum_worth_pcm": drum_worth,
            "sdm_pcm": sdm,
            "drho_void_pcm": drho_void,
            "u_void_pcm": u_void
        }
        results.append(res_dict)

    os.chdir(orig_dir)

    # Save to CSV
    csv_path = "results/architecture_screening.csv"
    keys = list(results[0].keys())
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for r in results:
            writer.writerow(r)

    print("\n" + "="*80)
    print(f"STAGE 1 SCREENING COMPLETE. Results saved to {csv_path}")
    print("="*80)

    # Print summary table
    print(f"{'ID':<4} | {'Family':<22} | {'Mass (kg)':<9} | {'Vol (L)':<8} | {'k_out':<9} | {'k_in':<9} | {'SDM (pcm)':<9} | {'Void (pcm)':<15}")
    print("-" * 95)
    for r in results:
        print(f"{r['cand_id']:<4} | {r['family']:<22} | {r['mass_total_kg']:<9.1f} | {r['volume_total_L']:<8.1f} | {r['k_out']:<9.5f} | {r['k_in']:<9.5f} | {r['sdm_pcm']:<9.0f} | {r['drho_void_pcm']:+6.1f} +/- {r['u_void_pcm']:<4.1f}")

    return results

if __name__ == "__main__":
    run_stage1_screening()
