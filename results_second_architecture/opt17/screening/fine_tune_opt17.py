#!/usr/bin/env python3
"""
FINE-TUNED SCREENING FOR OPT-17
Evaluate 4 targeted candidates aiming for strict no-regression across all metrics:
- Lifetime >= 28.0 years (k_op >= 1.0195)
- Mass <= 540.44 kg
- Drum Worth >= 4740.9 pcm
- SDM >= 2831.8 pcm
- SDM_stk >= 2393.7 pcm
- Fq < 3.428
- drho_void <= -87.1 pcm

Candidates:
1. Cand-FT1: 2-zone enrichment (18.2% inner 37 pins, 19.75% outer 180 pins), H_core=42.2 cm, ax_refl=9.40 cm, 150° drums, ax_duct_r=0.850 cm
2. Cand-FT2: 3-zone enrichment (18.0% rings 0-2 [19 pins], 19.0% rings 3-5 [72 pins], 19.75% rings 6-8 [126 pins]), H_core=42.3 cm, ax_refl=9.35 cm, 152° drums, ax_duct_r=0.860 cm
3. Cand-FT3: 2-zone enrichment (18.5% inner 37 pins, 19.75% outer 180 pins), pitch=2.270 cm, H_core=42.4 cm, ax_refl=9.30 cm, 150° drums, ax_duct_r=0.850 cm
4. Cand-FT4: 3-zone enrichment (18.2% inner 19 pins, 19.2% mid 72 pins, 19.75% outer 126 pins), H_core=42.2 cm, ax_refl=9.40 cm, 152° drums, ax_duct_r=0.860 cm
"""

import os
import sys
import math
import time
import csv
import numpy as np
import openmc

os.environ["PATH"] = "/home/aleri/miniforge3/envs/openmc-env/bin:" + os.environ.get("PATH", "")
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
os.environ["OMP_NUM_THREADS"] = "8"
openmc.config['cross_sections'] = "/home/aleri/wmp_data/cross_sections_wmp.xml"

SCREEN_DIR = "/home/aleri/OpenMC_projects/first_test/results_second_architecture/opt17/screening"
os.chdir(SCREEN_DIR)


class FTCandidateParams:
    def __init__(self, name,
                 pitch=2.270,
                 rf=0.580,
                 enr_zones=[(2, 18.0), (5, 19.0), (8, 19.75)], # (max_ring_in_zone, enrichment)
                 core_h=42.2, core_r=19.45,
                 refl_rad_th=9.80, refl_ax_th=9.40,
                 ax_duct_r=0.850,
                 drum_count=10, drum_r=3.95, abs_th=1.50, abs_angle=150.0,
                 r_add=0.150, r_add_pos=0.910):
        self.name = name
        self.pitch = pitch
        self.rf = rf
        self.enr_zones = enr_zones
        self.core_h = core_h
        self.core_r = core_r
        self.refl_rad_th = refl_rad_th
        self.refl_ax_th = refl_ax_th
        self.refl_rad_outer_r = core_r + refl_rad_th
        self.total_h = core_h + 2.0 * refl_ax_th
        self.ax_duct_r = ax_duct_r
        self.drum_count = drum_count
        self.drum_r = drum_r
        self.abs_th = abs_th
        self.abs_angle = abs_angle
        self.drum_clearance = 0.18
        self.drum_rad_pos = core_r + self.drum_clearance + drum_r
        self.r_add = r_add
        self.r_add_pos = r_add_pos
        self.clad_th = 0.050
        self.chan_gap = 0.075
        self.b10_enrichment = 0.90
        
        # Densities
        self.rho_un = 13.50
        self.rho_clad = 8.57
        self.rho_yh2 = 4.30
        self.rho_beo = 2.90
        self.rho_b4c = 2.52
        self.rho_he = 0.001146


def build_materials_ft(p: FTCandidateParams, is_void=False):
    openmc.reset_auto_ids()
    
    # Fuel materials for each zone
    fuel_mats = {}
    for zone_idx, (max_ring, enr) in enumerate(p.enr_zones):
        mat = openmc.Material(name=f"UN Fuel Zone {zone_idx} ({enr:.2f}%)")
        mat.set_density("g/cm3", p.rho_un)
        u235_at = (enr / 235.0) / (enr / 235.0 + (100.0 - enr) / 238.0)
        mat.add_nuclide("U235", 0.5 * u235_at, "ao")
        mat.add_nuclide("U238", 0.5 * (1.0 - u235_at), "ao")
        mat.add_nuclide("N14", 0.5 * 0.9963, "ao")
        mat.add_nuclide("N15", 0.5 * 0.0037, "ao")
        mat.temperature = 1050.0
        fuel_mats[zone_idx] = mat

    # Cladding
    clad = openmc.Material(name="Nb-1Zr Clad")
    clad.set_density("g/cm3", p.rho_clad)
    clad.add_element("Nb", 0.99, "wo")
    clad.add_element("Zr", 0.01, "wo")
    clad.temperature = 1050.0

    # Coolant
    cool = openmc.Material(name="Helium Coolant")
    cool_rho = 1.0e-7 if is_void else p.rho_he
    cool.set_density("g/cm3", cool_rho)
    cool.add_nuclide("He4", 1.0, "ao")
    cool.temperature = 1050.0

    # Moderator
    mod = openmc.Material(name="Canned YH2.0")
    mod.set_density("g/cm3", p.rho_yh2)
    mod.add_nuclide("Y89", 1.0, "ao")
    mod.add_nuclide("H1", 2.0, "ao")
    mod.add_s_alpha_beta("c_H_in_ZrH")
    mod.temperature = 800.0

    # Reflector
    refl = openmc.Material(name="BeO Reflector")
    refl.set_density("g/cm3", p.rho_beo)
    refl.add_nuclide("Be9", 1.0, "ao")
    refl.add_nuclide("O16", 1.0, "ao")
    refl.add_s_alpha_beta("c_Be_in_BeO")
    refl.add_s_alpha_beta("c_O_in_BeO")
    refl.temperature = 800.0

    # Absorber
    absorb = openmc.Material(name="B4C Absorber")
    absorb.set_density("g/cm3", p.rho_b4c)
    absorb.add_nuclide("B10", 4.0 * p.b10_enrichment, "ao")
    absorb.add_nuclide("B11", 4.0 * (1.0 - p.b10_enrichment), "ao")
    absorb.add_nuclide("C0", 1.0, "ao")
    absorb.temperature = 800.0

    all_mats = list(fuel_mats.values()) + [clad, cool, mod, refl, absorb]
    return openmc.Materials(all_mats), fuel_mats, clad, cool, mod, refl, absorb


def make_pin_univ_ft(p: FTCandidateParams, f_mat, clad, cool, mod, refl):
    rc_out = p.rf + p.clad_th
    rchan_out = rc_out + p.chan_gap

    z_core_bot = openmc.ZPlane(z0=-p.core_h/2.0)
    z_core_top = openmc.ZPlane(z0=+p.core_h/2.0)
    z_refl_bot = openmc.ZPlane(z0=-p.total_h/2.0, boundary_type="vacuum")
    z_refl_top = openmc.ZPlane(z0=+p.total_h/2.0, boundary_type="vacuum")

    cyl_f = openmc.ZCylinder(r=p.rf)
    cyl_c = openmc.ZCylinder(r=rc_out)
    cyl_ch = openmc.ZCylinder(r=rchan_out)

    aux_cyls = []
    for k in range(6):
        ang = math.radians(k * 60.0)
        x_c = p.r_add_pos * math.cos(ang)
        y_c = p.r_add_pos * math.sin(ang)
        aux_cyls.append(openmc.ZCylinder(x0=x_c, y0=y_c, r=p.r_add))

    c_fuel = openmc.Cell(fill=f_mat, region=-cyl_f & +z_core_bot & -z_core_top)
    c_clad = openmc.Cell(fill=clad, region=+cyl_f & -cyl_c & +z_core_bot & -z_core_top)
    c_cool_ann = openmc.Cell(fill=cool, region=+cyl_c & -cyl_ch & +z_core_bot & -z_core_top)

    aux_cells = []
    aux_union = None
    for c in aux_cyls:
        aux_cells.append(openmc.Cell(fill=cool, region=-c & +z_core_bot & -z_core_top))
        aux_union = -c if aux_union is None else (aux_union | -c)

    mod_reg = +cyl_ch & +z_core_bot & -z_core_top
    if aux_union is not None:
        mod_reg = mod_reg & ~aux_union
    c_mod = openmc.Cell(fill=mod, region=mod_reg)

    cyl_ax = openmc.ZCylinder(r=p.ax_duct_r)
    c_duct_bot = openmc.Cell(fill=cool, region=-cyl_ax & +z_refl_bot & -z_core_bot)
    c_duct_top = openmc.Cell(fill=cool, region=-cyl_ax & +z_core_top & -z_refl_top)

    aux_ducts_bot = []
    aux_ducts_top = []
    aux_ax_union = None
    for c in aux_cyls:
        aux_ducts_bot.append(openmc.Cell(fill=cool, region=-c & +z_refl_bot & -z_core_bot))
        aux_ducts_top.append(openmc.Cell(fill=cool, region=-c & +z_core_top & -z_refl_top))
        aux_ax_union = -c if aux_ax_union is None else (aux_ax_union | -c)

    c_beo_bot = openmc.Cell(fill=refl, region=+cyl_ax & +z_refl_bot & -z_core_bot & ~aux_ax_union)
    c_beo_top = openmc.Cell(fill=refl, region=+cyl_ax & +z_core_top & -z_refl_top & ~aux_ax_union)

    return openmc.Universe(cells=[c_fuel, c_clad, c_cool_ann, c_mod, c_duct_bot, c_beo_bot, c_duct_top, c_beo_top] + aux_cells + aux_ducts_bot + aux_ducts_top)


def build_geometry_ft(p: FTCandidateParams, fuel_mats, clad, cool, mod, refl, absorb, drum_angles=None):
    if drum_angles is None:
        drum_angles = [0.0] * p.drum_count

    # Build pin universe for each zone
    pin_universes = {}
    for z_idx, f_mat in fuel_mats.items():
        pin_universes[z_idx] = make_pin_univ_ft(p, f_mat, clad, cool, mod, refl)

    z_core_bot = openmc.ZPlane(z0=-p.core_h/2.0)
    z_core_top = openmc.ZPlane(z0=+p.core_h/2.0)
    z_refl_bot = openmc.ZPlane(z0=-p.total_h/2.0, boundary_type="vacuum")
    z_refl_top = openmc.ZPlane(z0=+p.total_h/2.0, boundary_type="vacuum")

    c_mod_buf = openmc.Cell(fill=mod, region=+z_core_bot & -z_core_top)
    c_beo_buf_bot = openmc.Cell(fill=refl, region=+z_refl_bot & -z_core_bot)
    c_beo_buf_top = openmc.Cell(fill=refl, region=+z_core_top & -z_refl_top)
    u_mod_buf = openmc.Universe(cells=[c_mod_buf, c_beo_buf_bot, c_beo_buf_top])

    num_rings = 9
    hex_lat = openmc.HexLattice()
    hex_lat.center = (0.0, 0.0)
    hex_lat.pitch = [p.pitch]
    hex_lat.orientation = 'y'

    rings = []
    pin_ring_map = {}
    pin_centers = []
    v1 = np.array([p.pitch, 0.0])
    v2 = np.array([p.pitch / 2.0, p.pitch * np.sqrt(3.0) / 2.0])

    for r in range(num_rings - 1, -1, -1):
        # Determine zone
        target_zone = 0
        for z_idx, (max_ring, _) in enumerate(p.enr_zones):
            if r <= max_ring:
                target_zone = z_idx
                break
        u_ring = pin_universes[target_zone]

        if r == 0:
            rings.append([u_ring])
            pin_centers.append((0.0, 0.0))
            pin_ring_map[(0.0, 0.0)] = 0
        else:
            ring_u = [u_ring] * (6 * r)
            rings.append(ring_u)
            corners = [r*v2, r*(v2-v1), -r*v1, -r*v2, -r*(v2-v1), r*v1]
            for corn in range(6):
                p_start = corners[corn]
                p_end = corners[(corn + 1) % 6]
                for step in range(r):
                    pt = p_start + (p_end - p_start) * (step / r)
                    pt_coords = (float(pt[0]), float(pt[1]))
                    pin_centers.append(pt_coords)
                    pin_ring_map[pt_coords] = r

    hex_lat.universes = rings
    hex_lat.outer = u_mod_buf

    cyl_core = openmc.ZCylinder(r=p.core_r)
    cyl_refl_outer = openmc.ZCylinder(r=p.refl_rad_outer_r, boundary_type="vacuum")

    cyl_drum_outer = openmc.ZCylinder(r=p.drum_r)
    cyl_drum_inner = openmc.ZCylinder(r=p.drum_r - p.abs_th)
    half_sec = np.radians(p.abs_angle / 2.0)
    plane_sec1 = openmc.Plane(a=np.sin(half_sec), b=-np.cos(half_sec), c=0.0, d=0.0)
    plane_sec2 = openmc.Plane(a=np.sin(half_sec), b=+np.cos(half_sec), c=0.0, d=0.0)
    plane_front = openmc.XPlane(x0=0.0)

    abs_region = -cyl_drum_outer & +cyl_drum_inner & +plane_sec1 & +plane_sec2 & +plane_front
    beo_drum_region = -cyl_drum_outer & ~abs_region
    u_drum = openmc.Universe(cells=[openmc.Cell(fill=absorb, region=abs_region),
                                    openmc.Cell(fill=refl, region=beo_drum_region)])

    drum_cells = []
    drum_regions = []
    for k in range(p.drum_count):
        phi_k = k * 2.0 * np.pi / p.drum_count
        x_k = p.drum_rad_pos * math.cos(phi_k)
        y_k = p.drum_rad_pos * math.sin(phi_k)
        cyl_k = openmc.ZCylinder(x0=x_k, y0=y_k, r=p.drum_r)
        
        rot_deg = k * (360.0 / p.drum_count) + drum_angles[k]
        c_dk = openmc.Cell(fill=u_drum, region=-cyl_k & +z_core_bot & -z_core_top)
        c_dk.translation = (x_k, y_k, 0.0)
        c_dk.rotation = (0.0, 0.0, rot_deg)
        drum_cells.append(c_dk)
        drum_regions.append(-cyl_k)

    drum_union = drum_regions[0]
    for reg in drum_regions[1:]:
        drum_union = drum_union | reg

    cell_core = openmc.Cell(fill=hex_lat, region=-cyl_core & +z_refl_bot & -z_refl_top)
    cell_refl_rad = openmc.Cell(fill=refl, region=+cyl_core & -cyl_refl_outer & +z_core_bot & -z_core_top & ~drum_union)
    cell_refl_bot = openmc.Cell(fill=refl, region=+cyl_core & -cyl_refl_outer & +z_refl_bot & -z_core_bot)
    cell_refl_top = openmc.Cell(fill=refl, region=+cyl_core & -cyl_refl_outer & +z_core_top & -z_refl_top)

    root = openmc.Universe(cells=[cell_core, cell_refl_rad, cell_refl_bot, cell_refl_top] + drum_cells)
    return openmc.Geometry(root), pin_centers, pin_ring_map


def compute_mass_ft(p: FTCandidateParams, pin_centers, pin_ring_map):
    n_pins = len(pin_centers)
    
    # Fuel mass and U-235
    v_f_pin = math.pi * p.rf**2 * p.core_h
    m_f_pin = v_f_pin * p.rho_un * 1e-3
    m_f_tot = n_pins * m_f_pin
    
    m_u235 = 0.0
    for c in pin_centers:
        r = pin_ring_map[c]
        enr = 19.75
        for max_ring, e in p.enr_zones:
            if r <= max_ring:
                enr = e
                break
        m_u_pin = m_f_pin * (238.0 / 252.0)
        m_u235 += m_u_pin * (enr / 100.0)

    # Clad
    v_c_pin = math.pi * ((p.rf + p.clad_th)**2 - p.rf**2) * p.core_h
    m_clad_kg = n_pins * v_c_pin * p.rho_clad * 1e-3

    # Coolant active
    v_ch_pin = math.pi * ((p.rf + p.clad_th + p.chan_gap)**2 - (p.rf + p.clad_th)**2) * p.core_h
    v_aux_core = 6 * math.pi * p.r_add**2 * p.core_h
    v_cool_core = n_pins * (v_ch_pin + v_aux_core)

    # Coolant axial
    v_duct_ax = math.pi * p.ax_duct_r**2 * (2.0 * p.refl_ax_th)
    v_aux_ax = 6 * math.pi * p.r_add**2 * (2.0 * p.refl_ax_th)
    v_cool_ax = n_pins * (v_duct_ax + v_aux_ax)

    # Moderator
    cell_area = (np.sqrt(3.0) / 2.0) * p.pitch**2
    v_cell_tot = n_pins * cell_area * p.core_h
    v_mod_tot = v_cell_tot - (n_pins * v_f_pin + n_pins * v_c_pin + v_cool_core)
    m_mod_kg = v_mod_tot * p.rho_yh2 * 1e-3

    # BeO
    v_beo_ax_gross = math.pi * p.refl_rad_outer_r**2 * (2.0 * p.refl_ax_th)
    v_beo_ax_net = v_beo_ax_gross - v_cool_ax
    m_beo_ax_kg = v_beo_ax_net * p.rho_beo * 1e-3

    v_cyl_active = math.pi * p.refl_rad_outer_r**2 * p.core_h
    v_refl_rad_gross = v_cyl_active - (math.pi * p.core_r**2 * p.core_h)
    v_all_drums = p.drum_count * (math.pi * p.drum_r**2 * p.core_h)
    v_refl_rad_stat = v_refl_rad_gross - v_all_drums

    v_abs_single = (p.abs_angle / 360.0) * math.pi * (p.drum_r**2 - (p.drum_r - p.abs_th)**2) * p.core_h
    v_abs_tot = p.drum_count * v_abs_single
    m_abs_kg = v_abs_tot * p.rho_b4c * 1e-3

    v_beo_drums = v_all_drums - v_abs_tot
    m_beo_rad_kg = (v_refl_rad_stat + v_beo_drums) * p.rho_beo * 1e-3
    m_beo_tot_kg = m_beo_ax_kg + m_beo_rad_kg

    m_tot = m_f_tot + m_clad_kg + m_mod_kg + m_beo_tot_kg + m_abs_kg
    v_tot_L = math.pi * p.refl_rad_outer_r**2 * p.total_h * 1e-3

    return {
        "m_total_kg": m_tot,
        "v_total_L": v_tot_L,
        "m_u235_kg": m_u235,
        "m_fuel_kg": m_f_tot,
        "m_clad_kg": m_clad_kg,
        "m_mod_kg": m_mod_kg,
        "m_beo_kg": m_beo_tot_kg,
        "m_abs_kg": m_abs_kg
    }


def eval_ft_candidate(cand_p: FTCandidateParams):
    print(f"\n=======================================================")
    print(f" EVALUATING FINE-TUNED CANDIDATE: {cand_p.name}")
    print(f"=======================================================")
    t0 = time.time()
    
    mats, fuel_mats, clad, cool, mod, refl, absorb = build_materials_ft(cand_p)
    geom_op, pin_centers, pin_map = build_geometry_ft(cand_p, fuel_mats, clad, cool, mod, refl, absorb, drum_angles=[0.0]*cand_p.drum_count)
    mass_info = compute_mass_ft(cand_p, pin_centers, pin_map)
    print(f"Mass: {mass_info['m_total_kg']:.2f} kg (Opt-16: 540.44 kg | Margin: {540.44 - mass_info['m_total_kg']:+.2f} kg)")
    print(f"Volume: {mass_info['v_total_L']:.2f} L | U-235: {mass_info['m_u235_kg']:.2f} kg")

    def make_settings(particles=50000, batches=40, inactive=10):
        s = openmc.Settings()
        s.run_mode = "eigenvalue"
        s.particles = particles
        s.batches = batches
        s.inactive = inactive
        s.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}
        space = openmc.stats.Box((-cand_p.core_r*0.6, -cand_p.core_r*0.6, -cand_p.core_h*0.3),
                                 (cand_p.core_r*0.6, cand_p.core_r*0.6, cand_p.core_h*0.3))
        s.source = openmc.IndependentSource(space=space)
        return s

    mesh = openmc.RegularMesh()
    mesh.dimension = [80, 80]
    mesh.lower_left = [-cand_p.refl_rad_outer_r, -cand_p.refl_rad_outer_r]
    mesh.upper_right = [+cand_p.refl_rad_outer_r, +cand_p.refl_rad_outer_r]
    mesh_filter = openmc.MeshFilter(mesh)
    t_fis = openmc.Tally(name="fission_mesh")
    t_fis.filters = [mesh_filter]
    t_fis.scores = ['fission']
    tallies = openmc.Tallies([t_fis])

    # 1. Operating State
    sett_op = make_settings(particles=60000, batches=45, inactive=10)
    model_op = openmc.Model(materials=mats, geometry=geom_op, settings=sett_op, tallies=tallies)
    model_op.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett_op.batches}.h5") as sp:
        k_op = float(sp.keff.nominal_value)
        s_op = float(sp.keff.std_dev)
        t_f = sp.get_tally(name="fission_mesh")
        fiss_data = t_f.get_values(scores=['fission']).reshape((80, 80))
    os.remove(f"statepoint.{sett_op.batches}.h5")

    # Peaking factor
    xs = np.linspace(-cand_p.refl_rad_outer_r, cand_p.refl_rad_outer_r, 80)
    ys = np.linspace(-cand_p.refl_rad_outer_r, cand_p.refl_rad_outer_r, 80)
    XX, YY = np.meshgrid(xs, ys)
    core_mask = (XX**2 + YY**2) <= (cand_p.core_r**2)
    fiss_core = fiss_data[core_mask]
    avg_fiss = np.mean(fiss_core[fiss_core > 0])
    max_fiss = np.max(fiss_core)
    Fq = max_fiss / avg_fiss if avg_fiss > 0 else 1.0
    excess_pcm = (k_op - 1.0) / k_op * 1e5
    burnup_rate = 68.2
    life_yr = excess_pcm / burnup_rate if excess_pcm > 0 else 0.0

    print(f"k_op: {k_op:.5f} +/- {s_op:.5f} | Excess: {excess_pcm:+.1f} pcm | Est Life: {life_yr:.1f} yr | Fq: {Fq:.3f} (Opt-16: 3.428)")

    # 2. Shutdown State (180° drums)
    geom_in, _, _ = build_geometry_ft(cand_p, fuel_mats, clad, cool, mod, refl, absorb, drum_angles=[180.0]*cand_p.drum_count)
    sett_in = make_settings(particles=40000, batches=35, inactive=10)
    model_in = openmc.Model(materials=mats, geometry=geom_in, settings=sett_in)
    model_in.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett_in.batches}.h5") as sp:
        k_in = float(sp.keff.nominal_value)
        s_in = float(sp.keff.std_dev)
    os.remove(f"statepoint.{sett_in.batches}.h5")
    worth_pcm = (1.0/k_in - 1.0/k_op) * 1e5
    sdm_pcm = (1.0 - k_in) / k_in * 1e5
    print(f"k_in: {k_in:.5f} +/- {s_in:.5f} | Worth: {worth_pcm:.1f} pcm (Opt-16: 4740.9) | SDM: {sdm_pcm:.1f} pcm (Opt-16: 2831.8)")

    # 3. Single Stuck Drum
    stk_angles = [0.0] + [180.0] * (cand_p.drum_count - 1)
    geom_stk, _, _ = build_geometry_ft(cand_p, fuel_mats, clad, cool, mod, refl, absorb, drum_angles=stk_angles)
    model_stk = openmc.Model(materials=mats, geometry=geom_stk, settings=sett_in)
    model_stk.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett_in.batches}.h5") as sp:
        k_stk = float(sp.keff.nominal_value)
        s_stk = float(sp.keff.std_dev)
    os.remove(f"statepoint.{sett_in.batches}.h5")
    sdm_stk_pcm = (1.0 - k_stk) / k_stk * 1e5
    print(f"k_stk: {k_stk:.5f} +/- {s_stk:.5f} | Stuck SDM: {sdm_stk_pcm:.1f} pcm (Opt-16: 2393.7)")

    # 4. Voided State (High-precision 100k particles x 45 batches)
    mats_vd, fuel_mats_vd, clad_vd, cool_vd, mod_vd, refl_vd, abs_vd = build_materials_ft(cand_p, is_void=True)
    geom_vd, _, _ = build_geometry_ft(cand_p, fuel_mats_vd, clad_vd, cool_vd, mod_vd, refl_vd, abs_vd, drum_angles=[0.0]*cand_p.drum_count)
    sett_vd = make_settings(particles=80000, batches=45, inactive=10)
    model_vd = openmc.Model(materials=mats_vd, geometry=geom_vd, settings=sett_vd)
    model_vd.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{sett_vd.batches}.h5") as sp:
        k_vd = float(sp.keff.nominal_value)
        s_vd = float(sp.keff.std_dev)
    os.remove(f"statepoint.{sett_vd.batches}.h5")

    drho_vd = (1.0/k_op - 1.0/k_vd) * 1e5
    u_drho_vd = math.sqrt((s_op/k_op**2)**2 + (s_vd/k_vd**2)**2) * 1e5
    b1_vd = drho_vd + 1.0 * u_drho_vd
    b2_vd = drho_vd + 2.0 * u_drho_vd
    b3_vd = drho_vd + 3.0 * u_drho_vd
    print(f"k_vd: {k_vd:.5f} +/- {s_vd:.5f} | Void drho: {drho_vd:+.1f} +/- {u_drho_vd:.1f} pcm (1s: {b1_vd:+.1f}, 2s: {b2_vd:+.1f}) (Opt-16: -87.1)")

    elapsed = time.time() - t0
    print(f"--> Evaluated in {elapsed:.1f} s")

    return {
        "candidate": cand_p.name,
        "mass_kg": mass_info['m_total_kg'],
        "volume_L": mass_info['v_total_L'],
        "k_op": k_op,
        "s_op": s_op,
        "excess_pcm": excess_pcm,
        "life_yr": life_yr,
        "k_in": k_in,
        "worth_pcm": worth_pcm,
        "sdm_pcm": sdm_pcm,
        "k_stk": k_stk,
        "sdm_stk_pcm": sdm_stk_pcm,
        "Fq": Fq,
        "k_vd": k_vd,
        "s_vd": s_vd,
        "drho_vd_pcm": drho_vd,
        "u_drho_vd_pcm": u_drho_vd,
        "b1_vd_pcm": b1_vd,
        "b2_vd_pcm": b2_vd,
        "b3_vd_pcm": b3_vd
    }


def main():
    candidates = [
        # Candidate FT1: 2-zone enrichment (18.5% inner 37 pins, 19.75% outer 180 pins), H_core=42.4 cm, ax_refl=9.30 cm, 150° drums, ax_duct_r=0.850 cm
        FTCandidateParams(
            name="Cand_FT1_2Zone_18.5pct",
            enr_zones=[(3, 18.50), (8, 19.75)],
            core_h=42.4, refl_ax_th=9.30,
            ax_duct_r=0.850, abs_angle=150.0, abs_th=1.50
        ),

        # Candidate FT2: 3-zone enrichment (18.0% 19 pins, 19.0% 72 pins, 19.75% 126 pins), H_core=42.5 cm, ax_refl=9.25 cm, 152° drums, ax_duct_r=0.860 cm
        FTCandidateParams(
            name="Cand_FT2_3Zone_18.0_19.0_19.75",
            enr_zones=[(2, 18.00), (5, 19.00), (8, 19.75)],
            core_h=42.5, refl_ax_th=9.25,
            ax_duct_r=0.860, abs_angle=152.0, abs_th=1.50
        ),

        # Candidate FT3: 2-zone enrichment (18.8% inner 37 pins, 19.75% outer 180 pins), H_core=42.2 cm, ax_refl=9.40 cm, 150° drums, ax_duct_r=0.850 cm
        FTCandidateParams(
            name="Cand_FT3_2Zone_18.8pct",
            enr_zones=[(3, 18.80), (8, 19.75)],
            core_h=42.2, refl_ax_th=9.40,
            ax_duct_r=0.850, abs_angle=150.0, abs_th=1.50
        ),

        # Candidate FT4: 3-zone enrichment (18.2% 19 pins, 19.2% 72 pins, 19.75% 126 pins), H_core=42.3 cm, ax_refl=9.35 cm, 152° drums, ax_duct_r=0.860 cm
        FTCandidateParams(
            name="Cand_FT4_3Zone_18.2_19.2_19.75",
            enr_zones=[(2, 18.20), (5, 19.20), (8, 19.75)],
            core_h=42.3, refl_ax_th=9.35,
            ax_duct_r=0.860, abs_angle=152.0, abs_th=1.50
        )
    ]

    results = []
    for c in candidates:
        res = eval_ft_candidate(c)
        results.append(res)

    csv_path = os.path.join(SCREEN_DIR, "fine_tuned_comparison.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0].keys()))
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved fine-tuned results to {csv_path}")
    print("\nSummary Comparison against Opt-16 Baseline:")
    print("Opt-16 Baseline: Mass=540.44kg, Life=28.0y, Worth=4740.9pcm, SDM=2831.8pcm, SDM_stk=2393.7pcm, Fq=3.428, Void=-87.1pcm")
    header = f"{'Candidate':<32} | {'Mass [kg]':<9} | {'Life':<6} | {'Worth':<8} | {'SDM [pcm]':<9} | {'SDM_stk':<9} | {'Fq':<6} | {'drho_vd [pcm]':<15}"
    print(header)
    print("-" * len(header))
    for r in results:
        vd_str = f"{r['drho_vd_pcm']:+.1f} +/- {r['u_drho_vd_pcm']:.1f}"
        print(f"{r['candidate']:<32} | {r['mass_kg']:9.2f} | {r['life_yr']:5.1f}y | {r['worth_pcm']:8.1f} | {r['sdm_pcm']:9.1f} | {r['sdm_stk_pcm']:9.1f} | {r['Fq']:6.3f} | {vd_str:<15}")


if __name__ == "__main__":
    main()
