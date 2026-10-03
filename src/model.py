"""
Conceptual HALEU-Fueled NEP Space Reactor Neutronics Model (Baseline 1 / Optimized)
===================================================================================
OpenMC Python API Model for Space Nuclear Electric Propulsion.

Features in Baseline 1 / Optimized:
1. Fuel: Sintered Uranium Nitride (UN) at 19.75 wt% U-235 (HALEU)
2. Cladding: Niobium-1% Zirconium (Nb-1Zr) refractory alloy
3. Coolant: Helium-Xenon (He-Xe, M=40 g/mol, 72% He, 28% Xe) in explicit annular channels
4. Moderator: Yttrium Hydride (YH2) matrix (acknowledging c_H_in_ZrH surrogate in ENDF/B-VII.1)
5. Reflectors: High-temperature Beryllium Oxide (BeO) radial and axial caps
6. Control: 6 external rotating BeO drums with 120-degree B4C absorber sectors
7. Full parameterization: Supports drum angle sweeps, stuck-drum safety, and coolant states.
"""

import os
import sys
import datetime
import numpy as np
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import openmc


# =============================================================================
# 1. PARAMETERS CONTAINER
# =============================================================================
class ReactorParameters:
    """Fully parameterized configuration container for Baseline 1 / Optimized NEP reactor."""

    def __init__(self,
                 fuel_radius=0.58,              # cm
                 cladding_thickness=0.05,       # cm (0.5 mm Nb-1Zr)
                 coolant_channel_thickness=0.045,# cm (0.45 mm annular channel)
                 moderator_fuel_pitch=1.615,    # cm (Opt-8 balanced under-moderated pitch)
                 core_radius=18.20,             # cm (Opt-8 reference radius)
                 active_core_height=41.80,      # cm (H/R = 2.30, H/D = 1.15)
                 radial_reflector_thickness=10.20,# cm (BeO radial reflector)
                 axial_reflector_thickness=9.50, # cm (BeO top & bottom blocks with streaming ducts)
                 number_of_control_drums=8,     # 8 rotating BeO control drums
                 control_drum_radius=4.40,      # cm
                 absorber_thickness=1.45,       # cm (enhanced B4C sector thickness)
                 absorber_sector_angular_width=140.0, # deg (enhanced B4C sector arc)
                 drum_clearance=0.20,           # clearance from core boundary [cm]
                 b10_enrichment=0.90,           # 90 at% B-10 enriched B4C for space control drums
                 fuel_enrichment=19.75,         # wt% U-235 (HALEU standard)
                 drum_angle_deg=0.0,            # 0 = drums out, 180 = drums in
                 coolant_state='operating',     # 'operating', 'cold', 'hot', 'voided'
                 particles=10000,
                 batches=35,
                 inactive=10):

        # Pin dimensions
        self.fuel_radius = fuel_radius
        self.cladding_thickness = cladding_thickness
        self.clad_outer_radius = self.fuel_radius + self.cladding_thickness # 0.63 cm
        self.coolant_channel_thickness = coolant_channel_thickness
        self.channel_outer_radius = self.clad_outer_radius + self.coolant_channel_thickness # 0.70 cm
        self.moderator_fuel_pitch = moderator_fuel_pitch

        # Core and reflector dimensions
        self.core_radius = core_radius
        self.active_core_height = active_core_height
        self.aspect_ratio_H_R = self.active_core_height / self.core_radius
        self.radial_reflector_thickness = radial_reflector_thickness
        self.axial_reflector_thickness = axial_reflector_thickness
        self.radial_reflector_outer_radius = self.core_radius + self.radial_reflector_thickness
        self.total_reactor_height = self.active_core_height + 2 * self.axial_reflector_thickness

        # Control drums
        self.number_of_control_drums = number_of_control_drums
        self.control_drum_radius = control_drum_radius
        self.absorber_thickness = absorber_thickness
        self.absorber_sector_angular_width = absorber_sector_angular_width
        self.drum_clearance = drum_clearance
        if self.drum_clearance is not None:
            self.drum_radial_position = self.core_radius + self.drum_clearance + self.control_drum_radius
            if self.drum_radial_position + self.control_drum_radius > self.radial_reflector_outer_radius:
                self.radial_reflector_outer_radius = self.drum_radial_position + self.control_drum_radius + 0.8
                self.radial_reflector_thickness = self.radial_reflector_outer_radius - self.core_radius
        else:
            self.drum_radial_position = (self.core_radius + self.radial_reflector_outer_radius) / 2.0
        self.drum_angle_deg = drum_angle_deg

        # Material densities & compositions
        self.fuel_enrichment = fuel_enrichment
        self.fuel_density = 13.50       # g/cm3
        self.cladding_density = 8.57   # g/cm3
        self.moderator_density = 4.30  # g/cm3 (YH2)
        self.reflector_density = 2.90  # g/cm3 (BeO)
        self.b10_enrichment = b10_enrichment
        self.absorber_density = 2.52   # g/cm3 (B4C)

        # Coolant thermodynamic state
        self.coolant_state = coolant_state
        if self.coolant_state == 'operating':
            self.coolant_density = 0.00916  # g/cm3 (2.0 MPa, 1050 K, M=40 g/mol)
            self.coolant_temp_K = 1050.0
        elif self.coolant_state == 'cold':
            self.coolant_density = 0.01069  # g/cm3 (2.0 MPa, 900 K)
            self.coolant_temp_K = 900.0
        elif self.coolant_state == 'hot':
            self.coolant_density = 0.00740  # g/cm3 (2.0 MPa, 1300 K)
            self.coolant_temp_K = 1300.0
        elif self.coolant_state == 'voided':
            self.coolant_density = 1.0e-7   # near zero (depressurized / vacuum)
            self.coolant_temp_K = 1050.0
        else:
            self.coolant_density = 0.00916
            self.coolant_temp_K = 1050.0

        # Monte Carlo settings
        self.particles = particles
        self.batches = batches
        self.inactive = inactive
        self.active = self.batches - self.inactive


# =============================================================================
# 2. MATERIALS DEFINITION
# =============================================================================
def build_materials(p: ReactorParameters):
    """Define all reactor materials including He-Xe coolant with explicit isotopes."""
    # 1. UN Fuel (19.75 wt% U-235 in U, natural N)
    fuel = openmc.Material(name="UN Fuel (19.75 wt% U-235)")
    fuel.set_density("g/cm3", p.fuel_density)
    fuel.add_nuclide("U235", 0.5 * 0.199522, "ao")
    fuel.add_nuclide("U238", 0.5 * 0.800478, "ao")
    fuel.add_nuclide("N14", 0.5 * 0.9963, "ao")
    fuel.add_nuclide("N15", 0.5 * 0.0037, "ao")

    # 2. Nb-1Zr Cladding (99 wt% Nb, 1 wt% Zr)
    cladding = openmc.Material(name="Nb-1Zr Cladding")
    cladding.set_density("g/cm3", p.cladding_density)
    cladding.add_element("Nb", 0.99, "wo")
    cladding.add_element("Zr", 0.01, "wo")

    # 3. He-Xe Binary Noble Gas Coolant (72 at% He, 28 at% Xe, M=40 g/mol)
    coolant = openmc.Material(name="He-Xe Coolant")
    coolant.set_density("g/cm3", p.coolant_density)
    coolant.add_nuclide("He4", 0.72, "ao")
    coolant.add_element("Xe", 0.28, "ao")

    # 4. YH2 Moderator (with c_H_in_ZrH surrogate)
    moderator = openmc.Material(name="YH2 Moderator")
    moderator.set_density("g/cm3", p.moderator_density)
    moderator.add_nuclide("Y89", 1.0, "ao")
    moderator.add_nuclide("H1", 2.0, "ao")
    moderator.add_s_alpha_beta("c_H_in_ZrH")

    # 5. BeO Reflector (S(alpha, beta) thermal scattering)
    reflector = openmc.Material(name="BeO Reflector")
    reflector.set_density("g/cm3", p.reflector_density)
    reflector.add_nuclide("Be9", 1.0, "ao")
    reflector.add_nuclide("O16", 1.0, "ao")
    reflector.add_s_alpha_beta("c_Be_in_BeO")
    reflector.add_s_alpha_beta("c_O_in_BeO")

    # 6. B4C Absorber (enriched or natural B)
    absorber = openmc.Material(name="B4C Absorber")
    absorber.set_density("g/cm3", p.absorber_density)
    if hasattr(p, 'b10_enrichment') and p.b10_enrichment is not None:
        absorber.add_nuclide("B10", 4.0 * p.b10_enrichment, "ao")
        absorber.add_nuclide("B11", 4.0 * (1.0 - p.b10_enrichment), "ao")
    else:
        absorber.add_element("B", 4.0, "ao")
    absorber.add_nuclide("C0", 1.0, "ao")

    materials = openmc.Materials([fuel, cladding, coolant, moderator, reflector, absorber])
    return materials, fuel, cladding, coolant, moderator, reflector, absorber


# =============================================================================
# 3. GEOMETRY CONSTRUCTION
# =============================================================================
def build_geometry(p: ReactorParameters, fuel, cladding, coolant, moderator, reflector, absorber, drum_angles=None):
    """
    Construct full 3-D OpenMC CSG geometry for Baseline 1 / Optimized NEP reactor.
    drum_angles: optional list of 6 angles [deg]. If None, all drums use p.drum_angle_deg.
    """
    if drum_angles is None:
        drum_angles = [p.drum_angle_deg] * p.number_of_control_drums

    # --- Pin Cell Universe with Annular He-Xe Coolant Passage ---
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

    # Lattice Grid Size
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
            # Check that entire pin assembly (including coolant channel) lies inside core radius
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

    # --- Control Drum Template Universe ---
    r_abs_in = p.control_drum_radius - p.absorber_thickness
    drum_outer_cyl = openmc.ZCylinder(r=p.control_drum_radius, name="Drum Outer Cyl")
    drum_inner_cyl = openmc.ZCylinder(r=r_abs_in, name="Drum Inner Cyl")

    # Analytical sector centered on +X axis
    half_phi_rad = np.radians(p.absorber_sector_angular_width / 2.0)
    sin_phi = np.sin(half_phi_rad)
    cos_phi = np.cos(half_phi_rad)
    p_sector_1 = openmc.Plane(a=sin_phi, b=-cos_phi, c=0.0, d=0.0)
    p_sector_2 = openmc.Plane(a=sin_phi, b=+cos_phi, c=0.0, d=0.0)
    p_x_pos = openmc.XPlane(x0=0.0)

    abs_region = -drum_outer_cyl & +drum_inner_cyl & +p_sector_1 & +p_sector_2 & +p_x_pos
    beo_region = -drum_outer_cyl & ~abs_region

    c_drum_abs = openmc.Cell(name="Drum B4C Sector", fill=absorber, region=abs_region)
    c_drum_beo = openmc.Cell(name="Drum BeO Body", fill=reflector, region=beo_region)
    u_drum = openmc.Universe(cells=[c_drum_abs, c_drum_beo], name="Control Drum Universe")

    # --- Global Surfaces & Bounds ---
    s_core = openmc.ZCylinder(r=p.core_radius, name="Core Outer Cylinder")
    s_outer = openmc.ZCylinder(r=p.radial_reflector_outer_radius, boundary_type="vacuum", name="Reactor Outer Cylinder")

    z_core_bot = openmc.ZPlane(z0=-p.active_core_height / 2, name="Core Bottom Plane")
    z_core_top = openmc.ZPlane(z0=p.active_core_height / 2, name="Core Top Plane")
    z_refl_bot = openmc.ZPlane(z0=-p.total_reactor_height / 2, boundary_type="vacuum", name="Reactor Bottom Plane")
    z_refl_top = openmc.ZPlane(z0=p.total_reactor_height / 2, boundary_type="vacuum", name="Reactor Top Plane")

    # --- 6 Symmetrically Positioned Control Drums ---
    drum_cells = []
    drum_regions = []
    drum_centers = []

    for k in range(p.number_of_control_drums):
        angle_rad = k * 2.0 * np.pi / p.number_of_control_drums
        drum_nominal_deg = k * 360.0 / p.number_of_control_drums
        xk = p.drum_radial_position * np.cos(angle_rad)
        yk = p.drum_radial_position * np.sin(angle_rad)
        drum_centers.append((xk, yk, drum_nominal_deg, drum_angles[k]))

        # Total orientation = drum nominal position angle + rotation angle
        total_rot_deg = drum_nominal_deg + drum_angles[k]

        dk = openmc.ZCylinder(x0=xk, y0=yk, r=p.control_drum_radius, name=f"Drum {k} Cylinder")
        cdk = openmc.Cell(name=f"Control Drum {k}", fill=u_drum,
                          region=-dk & +z_core_bot & -z_core_top)
        cdk.translation = (xk, yk, 0.0)
        cdk.rotation = (0.0, 0.0, total_rot_deg)

        drum_cells.append(cdk)
        drum_regions.append(-dk)

    drum_union = drum_regions[0]
    for reg in drum_regions[1:]:
        drum_union = drum_union | reg

    # --- Main Assembly Cells ---
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
    geometry = openmc.Geometry(root_universe)

    return geometry, intact_pin_centers, drum_centers


# =============================================================================
# 4. TALLIES DEFINITIONS
# =============================================================================
def build_tallies(p: ReactorParameters):
    """Define integral tallies and 2D spatial mesh tallies."""
    tallies = openmc.Tallies()

    # 1. Integral Physics Tally
    t_int = openmc.Tally(name="integral_physics")
    t_int.scores = ["flux", "fission", "absorption", "heating"]
    tallies.append(t_int)

    # 2. 2D Spatial Mesh Tally
    mesh = openmc.RegularMesh(name="radial_midplane_mesh")
    mesh.dimension = [70, 70, 1]
    mesh.lower_left = [-p.radial_reflector_outer_radius, -p.radial_reflector_outer_radius, -p.active_core_height / 2]
    mesh.upper_right = [p.radial_reflector_outer_radius, p.radial_reflector_outer_radius, p.active_core_height / 2]

    mf = openmc.MeshFilter(mesh)
    t_mesh = openmc.Tally(name="spatial_mesh_tally")
    t_mesh.filters = [mf]
    t_mesh.scores = ["flux", "fission"]
    tallies.append(t_mesh)

    return tallies, mesh


# =============================================================================
# 5. SIMULATION SETTINGS
# =============================================================================
def build_settings(p: ReactorParameters):
    """Define Monte Carlo eigenvalue calculation settings."""
    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = p.particles
    settings.batches = p.batches
    settings.inactive = p.inactive

    box_w = p.core_radius * 0.65
    box_z = p.active_core_height * 0.35
    source_space = openmc.stats.Box(
        lower_left=[-box_w, -box_w, -box_z],
        upper_right=[box_w, box_w, box_z]
    )
    settings.source = openmc.IndependentSource(space=source_space)
    settings.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}
    return settings


# =============================================================================
# 6. EXACT ANALYTICAL MASS AND VOLUME INVENTORY
# =============================================================================
def calculate_inventory(p: ReactorParameters, num_pins: int):
    """
    Rigorously calculate analytical mass and volume inventory of all modeled
    components based on exact CSG geometry parameters.
    """
    # 1. UN Fuel
    v_fuel_per_pin = np.pi * (p.fuel_radius**2) * p.active_core_height
    v_fuel_total = num_pins * v_fuel_per_pin
    m_fuel_total = v_fuel_total * p.fuel_density

    # 2. Nb-1Zr Cladding
    v_clad_per_pin = np.pi * (p.clad_outer_radius**2 - p.fuel_radius**2) * p.active_core_height
    v_clad_total = num_pins * v_clad_per_pin
    m_clad_total = v_clad_total * p.cladding_density

    # 3. He-Xe Coolant in Annular Channels
    v_cool_per_pin = np.pi * (p.channel_outer_radius**2 - p.clad_outer_radius**2) * p.active_core_height
    v_cool_total = num_pins * v_cool_per_pin
    m_cool_total = v_cool_total * p.coolant_density

    # 4. YH2 Moderator Matrix (Active core volume minus pin channel assemblies)
    v_active_core = np.pi * (p.core_radius**2) * p.active_core_height
    v_pin_assemblies = num_pins * np.pi * (p.channel_outer_radius**2) * p.active_core_height
    v_mod_total = v_active_core - v_pin_assemblies
    m_mod_total = v_mod_total * p.moderator_density

    # 5. Control Drums (6 drums)
    r_abs_in = p.control_drum_radius - p.absorber_thickness
    area_abs_sector = (p.absorber_sector_angular_width / 360.0) * np.pi * (p.control_drum_radius**2 - r_abs_in**2)
    v_abs_1_drum = area_abs_sector * p.active_core_height
    v_abs_total = p.number_of_control_drums * v_abs_1_drum
    m_abs_total = v_abs_total * p.absorber_density

    v_drum_cyl_1 = np.pi * (p.control_drum_radius**2) * p.active_core_height
    v_drum_beo_1 = v_drum_cyl_1 - v_abs_1_drum
    v_drum_beo_total = p.number_of_control_drums * v_drum_beo_1
    m_drum_beo_total = v_drum_beo_total * p.reflector_density

    # 6. Radial Reflector Matrix (Annulus minus 6 drum cylinders)
    area_rad_annulus = np.pi * (p.radial_reflector_outer_radius**2 - p.core_radius**2)
    v_rad_annulus = area_rad_annulus * p.active_core_height
    v_rad_matrix_beo = v_rad_annulus - (p.number_of_control_drums * v_drum_cyl_1)
    m_rad_matrix_beo = v_rad_matrix_beo * p.reflector_density

    v_rad_refl_total = v_drum_beo_total + v_rad_matrix_beo
    m_rad_refl_total = v_rad_refl_total * p.reflector_density

    # 7. Axial Reflectors (Top and Bottom full-radius BeO caps)
    area_full_reactor = np.pi * (p.radial_reflector_outer_radius**2)
    v_axial_refl_total = 2.0 * area_full_reactor * p.axial_reflector_thickness
    m_axial_refl_total = v_axial_refl_total * p.reflector_density

    # Totals
    total_volume_cm3 = np.pi * (p.radial_reflector_outer_radius**2) * p.total_reactor_height
    total_mass_g = (m_fuel_total + m_clad_total + m_cool_total + m_mod_total +
                    m_rad_matrix_beo + m_drum_beo_total + m_axial_refl_total + m_abs_total)
    total_mass_kg = total_mass_g / 1000.0

    items = [
        {"name": "UN Fuel (19.75 wt% U-235)", "density": p.fuel_density,
         "volume_cm3": v_fuel_total, "volume_L": v_fuel_total / 1000.0,
         "mass_g": m_fuel_total, "mass_kg": m_fuel_total / 1000.0},
        {"name": "Nb-1Zr Fuel Cladding", "density": p.cladding_density,
         "volume_cm3": v_clad_total, "volume_L": v_clad_total / 1000.0,
         "mass_g": m_clad_total, "mass_kg": m_clad_total / 1000.0},
        {"name": "He-Xe Coolant (Channels)", "density": p.coolant_density,
         "volume_cm3": v_cool_total, "volume_L": v_cool_total / 1000.0,
         "mass_g": m_cool_total, "mass_kg": m_cool_total / 1000.0},
        {"name": "YH2 Moderator", "density": p.moderator_density,
         "volume_cm3": v_mod_total, "volume_L": v_mod_total / 1000.0,
         "mass_g": m_mod_total, "mass_kg": m_mod_total / 1000.0},
        {"name": "BeO Radial Reflector", "density": p.reflector_density,
         "volume_cm3": v_rad_refl_total, "volume_L": v_rad_refl_total / 1000.0,
         "mass_g": m_rad_refl_total, "mass_kg": m_rad_refl_total / 1000.0},
        {"name": "BeO Axial Reflectors (Top+Bot)", "density": p.reflector_density,
         "volume_cm3": v_axial_refl_total, "volume_L": v_axial_refl_total / 1000.0,
         "mass_g": m_axial_refl_total, "mass_kg": m_axial_refl_total / 1000.0},
        {"name": "B4C Absorber (Control Drums)", "density": p.absorber_density,
         "volume_cm3": v_abs_total, "volume_L": v_abs_total / 1000.0,
         "mass_g": m_abs_total, "mass_kg": m_abs_total / 1000.0},
    ]

    for it in items:
        it["mass_pct"] = (it["mass_g"] / total_mass_g) * 100.0
        it["vol_pct"] = (it["volume_cm3"] / total_volume_cm3) * 100.0

    summary = {
        "items": items,
        "total_volume_cm3": total_volume_cm3,
        "total_volume_L": total_volume_cm3 / 1000.0,
        "total_mass_g": total_mass_g,
        "total_mass_kg": total_mass_kg,
        "total_fuel_volume_L": v_fuel_total / 1000.0,
        "total_fuel_mass_kg": m_fuel_total / 1000.0,
        "total_moderator_volume_L": v_mod_total / 1000.0,
        "total_moderator_mass_kg": m_mod_total / 1000.0,
        "total_reflector_volume_L": (v_rad_refl_total + v_axial_refl_total) / 1000.0,
        "total_reflector_mass_kg": (m_rad_refl_total + m_axial_refl_total) / 1000.0,
        "total_absorber_mass_kg": m_abs_total / 1000.0,
        "num_pins": num_pins
    }
    return summary


# =============================================================================
# 7. REQUIRED VISUALIZATIONS GENERATOR
# =============================================================================
def generate_all_plots(p: ReactorParameters, intact_pin_centers, drum_centers, flux_slice, fiss_slice, output_dir):
    """
    Generate all required high-resolution plots:
    1. results/reactor_geometry.png (radial midplane view)
    2. results/reactor_vertical.png (vertical longitudinal cross-section)
    3. results/fuel_cell_zoom.png (close-up of 3 adjacent fuel pins with channels and clad)
    4. results/flux_distribution.png (spatial flux heatmap)
    5. results/fission_rate_distribution.png (spatial fission heatmap)
    """
    os.makedirs(output_dir, exist_ok=True)

    # -------------------------------------------------------------------------
    # Plot 1: Radial Midplane Geometry
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 8), dpi=300)
    ax.set_aspect('equal')
    ax.set_xlim(-p.radial_reflector_outer_radius - 2, p.radial_reflector_outer_radius + 2)
    ax.set_ylim(-p.radial_reflector_outer_radius - 2, p.radial_reflector_outer_radius + 2)
    ax.set_xlabel("X Position [cm]", fontsize=12)
    ax.set_ylabel("Y Position [cm]", fontsize=12)
    ax.set_title(f"HALEU NEP Optimized Reactor - Radial Midplane (H/R={p.aspect_ratio_H_R:.2f})", fontsize=12, fontweight='bold')

    # Radial Reflector Annulus (BeO)
    refl_patch = patches.Circle((0, 0), p.radial_reflector_outer_radius,
                                facecolor='#C8E6C9', edgecolor='#2E7D32', linewidth=2, label="BeO Reflector")
    ax.add_patch(refl_patch)

    # Active Core Region (YH2 Moderator)
    core_patch = patches.Circle((0, 0), p.core_radius,
                                facecolor='#BBDEFB', edgecolor='#1565C0', linewidth=1.5, label="YH2 Moderator")
    ax.add_patch(core_patch)

    # Fuel Pins with Cladding and He-Xe Coolant
    for x, y in intact_pin_centers:
        ch_p = patches.Circle((x, y), p.channel_outer_radius, facecolor='#FFE082', edgecolor='none')
        cl_p = patches.Circle((x, y), p.clad_outer_radius, facecolor='#B0BEC5', edgecolor='none')
        f_p = patches.Circle((x, y), p.fuel_radius, facecolor='#D32F2F', edgecolor='none')
        ax.add_patch(ch_p)
        ax.add_patch(cl_p)
        ax.add_patch(f_p)

    ax.scatter([], [], c='#D32F2F', s=50, label=f"UN Fuel ({len(intact_pin_centers)} Pins)")
    ax.scatter([], [], c='#B0BEC5', s=70, label="Nb-1Zr Clad (0.5 mm)")
    ax.scatter([], [], c='#FFE082', s=90, label="He-Xe Coolant Channel (0.7 mm)")

    # Control Drums
    for k, (xd, yd, nominal_deg, rot_deg) in enumerate(drum_centers):
        d_beo = patches.Circle((xd, yd), p.control_drum_radius, facecolor='#81C784', edgecolor='#1B5E20', linewidth=1.2)
        ax.add_patch(d_beo)
        total_deg = nominal_deg + rot_deg
        half_arc = p.absorber_sector_angular_width / 2.0
        t1 = total_deg - half_arc
        t2 = total_deg + half_arc
        d_abs = patches.Wedge((xd, yd), p.control_drum_radius, t1, t2,
                              width=p.absorber_thickness, facecolor='#212121', edgecolor='#000000', linewidth=1.2)
        ax.add_patch(d_abs)

    ax.scatter([], [], c='#212121', marker='s', s=80, label=f"B4C Absorber ({p.absorber_sector_angular_width:.0f}° Sector)")
    ax.scatter([], [], c='#81C784', s=100, label="BeO Control Drum")

    core_bound = patches.Circle((0, 0), p.core_radius, facecolor='none', edgecolor='#0D47A1', linestyle='--', linewidth=1.5)
    ax.add_patch(core_bound)

    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper right', framealpha=0.95, fontsize=8.5)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "reactor_geometry.png"))
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Plot 2: Full Vertical Longitudinal Cross-Section
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 9), dpi=300)
    ax.set_aspect('equal')
    ax.set_xlim(-p.radial_reflector_outer_radius - 6, p.radial_reflector_outer_radius + 6)
    ax.set_ylim(-p.total_reactor_height / 2 - 5, p.total_reactor_height / 2 + 5)
    ax.set_xlabel("Radial Distance from Core Axis X [cm]", fontsize=12)
    ax.set_ylabel("Axial Position Z [cm]", fontsize=12)
    ax.set_title("HALEU NEP Reactor - Full Vertical Longitudinal Section", fontsize=12, fontweight='bold')

    # Top and Bottom Axial Reflectors
    z_bot = -p.total_reactor_height / 2
    z_core_bot = -p.active_core_height / 2
    z_core_top = p.active_core_height / 2
    z_top = p.total_reactor_height / 2
    r_refl = p.radial_reflector_outer_radius
    r_core = p.core_radius

    # Bottom BeO reflector
    bot_refl_patch = patches.Rectangle((-r_refl, z_bot), 2 * r_refl, p.axial_reflector_thickness,
                                       facecolor='#C8E6C9', edgecolor='#2E7D32', linewidth=1.5, label="BeO Axial Reflector")
    ax.add_patch(bot_refl_patch)

    # Top BeO reflector
    top_refl_patch = patches.Rectangle((-r_refl, z_core_top), 2 * r_refl, p.axial_reflector_thickness,
                                       facecolor='#C8E6C9', edgecolor='#2E7D32', linewidth=1.5)
    ax.add_patch(top_refl_patch)

    # Radial Reflector sides
    left_refl = patches.Rectangle((-r_refl, z_core_bot), p.radial_reflector_thickness, p.active_core_height,
                                  facecolor='#E8F5E9', edgecolor='#2E7D32', linewidth=1.5, label="BeO Radial Reflector")
    right_refl = patches.Rectangle((r_core, z_core_bot), p.radial_reflector_thickness, p.active_core_height,
                                   facecolor='#E8F5E9', edgecolor='#2E7D32', linewidth=1.5)
    ax.add_patch(left_refl)
    ax.add_patch(right_refl)

    # Active Core Matrix (YH2 Moderator)
    core_vert = patches.Rectangle((-r_core, z_core_bot), 2 * r_core, p.active_core_height,
                                  facecolor='#BBDEFB', edgecolor='#1565C0', linewidth=1.5, label="YH2 Moderator Matrix")
    ax.add_patch(core_vert)

    # Vertical Fuel Pins in Slice (showing pins near y=0)
    for x, y in intact_pin_centers:
        if abs(y) < p.moderator_fuel_pitch / 2:
            # Channel
            ch_v = patches.Rectangle((x - p.channel_outer_radius, z_core_bot), 2 * p.channel_outer_radius, p.active_core_height,
                                     facecolor='#FFE082', edgecolor='none')
            cl_v = patches.Rectangle((x - p.clad_outer_radius, z_core_bot), 2 * p.clad_outer_radius, p.active_core_height,
                                     facecolor='#B0BEC5', edgecolor='none')
            f_v = patches.Rectangle((x - p.fuel_radius, z_core_bot), 2 * p.fuel_radius, p.active_core_height,
                                    facecolor='#D32F2F', edgecolor='none')
            ax.add_patch(ch_v)
            ax.add_patch(cl_v)
            ax.add_patch(f_v)

    # Control Drums in cross-section
    for k, (xd, yd, nominal_deg, rot_deg) in enumerate(drum_centers):
        if abs(yd) < p.control_drum_radius:
            d_rect = patches.Rectangle((xd - p.control_drum_radius, z_core_bot), 2 * p.control_drum_radius, p.active_core_height,
                                       facecolor='#81C784', edgecolor='#1B5E20', linewidth=1.2)
            ax.add_patch(d_rect)
            # Absorber strip
            if xd > 0:
                abs_rect = patches.Rectangle((xd + p.control_drum_radius - p.absorber_thickness, z_core_bot),
                                             p.absorber_thickness, p.active_core_height, facecolor='#212121')
            else:
                abs_rect = patches.Rectangle((xd - p.control_drum_radius, z_core_bot),
                                             p.absorber_thickness, p.active_core_height, facecolor='#212121')
            ax.add_patch(abs_rect)

    # Dimension arrows & annotations
    # Core height
    ax.annotate('', xy=(-r_refl - 2, z_core_top), xytext=(-r_refl - 2, z_core_bot),
                arrowprops=dict(arrowstyle='<->', color='navy', lw=1.5))
    ax.text(-r_refl - 4, 0, f"Active Height\nH = {p.active_core_height:.1f} cm", rotation=90, va='center', ha='center', fontsize=9, color='navy')

    # Total height
    ax.annotate('', xy=(r_refl + 3, z_top), xytext=(r_refl + 3, z_bot),
                arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))
    ax.text(r_refl + 5, 0, f"Total Height\nH_tot = {p.total_reactor_height:.1f} cm", rotation=90, va='center', ha='center', fontsize=9)

    # Core Diameter
    ax.annotate('', xy=(-r_core, z_core_bot - 1.5), xytext=(r_core, z_core_bot - 1.5),
                arrowprops=dict(arrowstyle='<->', color='blue', lw=1.5))
    ax.text(0, z_core_bot - 3.5, f"Core Diameter D = {2*p.core_radius:.1f} cm (R={p.core_radius:.1f} cm)", ha='center', fontsize=9, color='blue')

    # Reflector thickness
    ax.annotate('', xy=(r_core, z_core_top + 1.5), xytext=(r_refl, z_core_top + 1.5),
                arrowprops=dict(arrowstyle='<->', color='darkgreen', lw=1.5))
    ax.text((r_core + r_refl)/2, z_core_top + 3.0, f"t_refl = {p.radial_reflector_thickness:.1f} cm", ha='center', fontsize=8.5, color='darkgreen')

    ax.scatter([], [], c='#D32F2F', s=50, label="UN Fuel Pin")
    ax.scatter([], [], c='#FFE082', s=70, label="He-Xe Gas Flow Channel")
    ax.scatter([], [], c='#81C784', s=90, label="Control Drum (BeO + B4C)")

    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper right', framealpha=0.95, fontsize=8.5)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "reactor_vertical.png"))
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Plot 3: Fuel Cell Microscopic Close-Up
    # -------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 7), dpi=300)
    ax.set_aspect('equal')
    p_pitch = p.moderator_fuel_pitch
    ax.set_xlim(-1.7 * p_pitch, 1.7 * p_pitch)
    ax.set_ylim(-1.1 * p_pitch, 1.1 * p_pitch)
    ax.set_xlabel("X Position [cm]", fontsize=12)
    ax.set_ylabel("Y Position [cm]", fontsize=12)
    ax.set_title("HALEU NEP Fuel Cell Architecture - Microscopic Zoom", fontsize=12, fontweight='bold')

    # Moderator background
    bg_box = patches.Rectangle((-1.7 * p_pitch, -1.1 * p_pitch), 3.4 * p_pitch, 2.2 * p_pitch,
                               facecolor='#BBDEFB', edgecolor='#1565C0', linewidth=1.5, label="YH2 Moderator Matrix")
    ax.add_patch(bg_box)

    # 3 adjacent pins at (-p, 0), (0, 0), (+p, 0)
    pin_xs = [-p_pitch, 0.0, p_pitch]
    for px in pin_xs:
        # Coolant channel
        c_ch = patches.Circle((px, 0), p.channel_outer_radius, facecolor='#FFE082', edgecolor='#F57F17', linewidth=1.2)
        # Cladding
        c_cl = patches.Circle((px, 0), p.clad_outer_radius, facecolor='#B0BEC5', edgecolor='#37474F', linewidth=1.2)
        # Fuel pellet
        c_fu = patches.Circle((px, 0), p.fuel_radius, facecolor='#D32F2F', edgecolor='#B71C1C', linewidth=1.2)
        ax.add_patch(c_ch)
        ax.add_patch(c_cl)
        ax.add_patch(c_fu)

    # Lattice unit cell boundary around center pin
    unit_box = patches.Rectangle((-p_pitch/2, -p_pitch/2), p_pitch, p_pitch,
                                 facecolor='none', edgecolor='#0D47A1', linestyle='--', linewidth=1.8, label=f"Lattice Unit Cell (Pitch={p_pitch:.2f} cm)")
    ax.add_patch(unit_box)

    # Callout Annotations on central pin
    ax.annotate(f"UN Fuel Pellet\nr_fuel = {p.fuel_radius*10:.1f} mm\n(19.75 wt% U-235)",
                xy=(0, 0.2), xytext=(-0.8 * p_pitch, 0.7 * p_pitch),
                arrowprops=dict(arrowstyle='->', color='#B71C1C', lw=1.5), fontsize=9, fontweight='bold', color='#B71C1C')

    ax.annotate(f"Nb-1Zr Cladding\nt_clad = {p.cladding_thickness*10:.1f} mm\n(r_clad = {p.clad_outer_radius*10:.1f} mm)",
                xy=(0.60, 0.2), xytext=(0.4 * p_pitch, 0.75 * p_pitch),
                arrowprops=dict(arrowstyle='->', color='#37474F', lw=1.5), fontsize=9, color='#37474F')

    ax.annotate(f"He-Xe Coolant Gap\nt_gap = {p.coolant_channel_thickness*10:.1f} mm (2.0 MPa, 1050 K)\n(r_channel = {p.channel_outer_radius*10:.1f} mm)",
                xy=(0.66, -0.2), xytext=(0.3 * p_pitch, -0.85 * p_pitch),
                arrowprops=dict(arrowstyle='->', color='#F57F17', lw=1.5), fontsize=9, color='#E65100')

    # Pitch annotation
    ax.annotate('', xy=(-p_pitch, -0.8 * p_pitch), xytext=(0, -0.8 * p_pitch),
                arrowprops=dict(arrowstyle='<->', color='navy', lw=1.5))
    ax.text(-p_pitch/2, -0.98 * p_pitch, f"Pitch = {p_pitch*10:.1f} mm", ha='center', fontsize=9, color='navy')

    ax.scatter([], [], c='#D32F2F', s=60, label="UN Fuel Pellet")
    ax.scatter([], [], c='#B0BEC5', s=80, label="Nb-1Zr Cladding")
    ax.scatter([], [], c='#FFE082', s=100, label="He-Xe Coolant Channel")

    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(loc='upper right', framealpha=0.95, fontsize=8.5)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "fuel_cell_zoom.png"))
    plt.close(fig)

    # -------------------------------------------------------------------------
    # Plot 4 & 5: Flux and Fission Maps
    # -------------------------------------------------------------------------
    extent = [-p.radial_reflector_outer_radius, p.radial_reflector_outer_radius,
              -p.radial_reflector_outer_radius, p.radial_reflector_outer_radius]

    # Flux
    fig, ax = plt.subplots(figsize=(8, 7), dpi=300)
    im = ax.imshow(flux_slice, origin='lower', extent=extent, cmap='viridis')
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Relative Neutron Flux [particles/cm²·source]", fontsize=11)
    ax.add_patch(patches.Circle((0, 0), p.core_radius, facecolor='none', edgecolor='red', linestyle='--', linewidth=1.8, label=f"Active Core (R={p.core_radius} cm)"))
    ax.add_patch(patches.Circle((0, 0), p.radial_reflector_outer_radius, facecolor='none', edgecolor='white', linestyle='-', linewidth=1.5, label=f"Reflector Outer (R={p.radial_reflector_outer_radius} cm)"))
    ax.set_xlabel("X Position [cm]", fontsize=12)
    ax.set_ylabel("Y Position [cm]", fontsize=12)
    ax.set_title("HALEU NEP Core - Spatial Neutron Flux Distribution", fontsize=12, fontweight='bold')
    ax.legend(loc='upper right', framealpha=0.9, fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "flux_distribution.png"))
    plt.close(fig)

    # Fission
    fig, ax = plt.subplots(figsize=(8, 7), dpi=300)
    im = ax.imshow(fiss_slice, origin='lower', extent=extent, cmap='inferno')
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Relative Fission Rate [fissions/source]", fontsize=11)
    ax.add_patch(patches.Circle((0, 0), p.core_radius, facecolor='none', edgecolor='cyan', linestyle='--', linewidth=1.8, label=f"Active Core (R={p.core_radius} cm)"))
    ax.set_xlabel("X Position [cm]", fontsize=12)
    ax.set_ylabel("Y Position [cm]", fontsize=12)
    ax.set_title("HALEU NEP Core - Spatial Fission Rate Distribution", fontsize=12, fontweight='bold')
    ax.legend(loc='upper right', framealpha=0.9, fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, "fission_rate_distribution.png"))
    plt.close(fig)

    print("Generated all core visualization plots successfully!")


# =============================================================================
# 8. MAIN STANDALONE EXECUTION
# =============================================================================
def main():
    p = ReactorParameters()
    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 75)
    print("HALEU NEP REACTOR MODEL (BASELINE 1 / OPTIMIZED)")
    print("=" * 75)

    mats, fuel, clad, cool, mod, beo, b4c = build_materials(p)
    geom, intact_pins, drum_centers = build_geometry(p, fuel, clad, cool, mod, beo, b4c)
    tallies, mesh = build_tallies(p)
    settings = build_settings(p)

    model = openmc.Model(materials=mats, geometry=geom, settings=settings, tallies=tallies)
    model.export_to_xml()
    print("Exported XML files. Running simulation...")
    openmc.run()

    sp_file = f"statepoint.{p.batches}.h5"
    with openmc.StatePoint(sp_file) as sp:
        keff = sp.keff
        mesh_t = sp.get_tally(name="spatial_mesh_tally")
        dim = mesh.dimension
        flux_slice = mesh_t.get_slice(scores=["flux"]).mean.reshape((dim[0], dim[1]))
        fiss_slice = mesh_t.get_slice(scores=["fission"]).mean.reshape((dim[0], dim[1]))

    print(f"Simulation completed! keff = {keff.nominal_value:.5f} +/- {keff.std_dev:.5f}")

    inv = calculate_inventory(p, len(intact_pins))
    generate_all_plots(p, intact_pins, drum_centers, flux_slice, fiss_slice, output_dir)

    # Export material_inventory.csv
    with open(os.path.join(output_dir, "material_inventory.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Material Name", "Density [g/cm3]", "Volume [cm3]", "Volume [L]", "Mass [g]", "Mass [kg]", "Mass Fraction [%]", "Volume Fraction [%]"])
        for it in inv["items"]:
            writer.writerow([it["name"], f"{it['density']:.2f}", f"{it['volume_cm3']:.2f}", f"{it['volume_L']:.4f}", f"{it['mass_g']:.2f}", f"{it['mass_kg']:.3f}", f"{it['mass_pct']:.2f}", f"{it['vol_pct']:.2f}"])
        writer.writerow([])
        writer.writerow(["Total Modeled Reactor", "-", f"{inv['total_volume_cm3']:.2f}", f"{inv['total_volume_L']:.4f}", f"{inv['total_mass_g']:.2f}", f"{inv['total_mass_kg']:.3f}", "100.00", "100.00"])

    # Export tallies.csv
    with open(os.path.join(output_dir, "tallies.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Tally Score", "Mean [per source particle]", "Std Dev", "Relative Error [%]"])
        with openmc.StatePoint(sp_file) as sp:
            t_int = sp.get_tally(name="integral_physics")
            for score in t_int.scores:
                val = t_int.get_slice(scores=[score]).mean[0, 0, 0]
                std = t_int.get_slice(scores=[score]).std_dev[0, 0, 0]
                rel = (std / val * 100.0) if val != 0 else 0.0
                writer.writerow([score, f"{val:.6e}", f"{std:.6e}", f"{rel:.2f}"])

    return p, inv, keff


if __name__ == "__main__":
    main()
