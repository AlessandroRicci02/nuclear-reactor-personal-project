#!/usr/bin/env python3
"""
TEST CANDIDATE OPT-17 CHAMPION CANDIDATE
Tests:
- Mass inventory & volume
- k_op (0° drums OUT), excess reactivity, estimated lifetime
- k_in (180° drums IN), drum worth, SDM
- k_stk (drum 0 OUT, drums 1-9 IN), single stuck drum SDM
- Fq radial power peaking factor
- drho_void (vacuum depressurization)
"""

import os
import sys
import math
import time
import numpy as np
import openmc

OPENMC_EXEC = "/home/aleri/miniforge3/envs/openmc-env/bin/openmc"
os.environ["PATH"] = "/home/aleri/miniforge3/envs/openmc-env/bin:" + os.environ.get("PATH", "")
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
os.environ["OMP_NUM_THREADS"] = "8"
openmc.config['cross_sections'] = "/home/aleri/wmp_data/cross_sections_wmp.xml"

SCREEN_DIR = "/home/aleri/OpenMC_projects/first_test/results_second_architecture/opt17/screening"
os.chdir(SCREEN_DIR)

# Core & Geometry parameters
pitch = 2.270
rf = 0.580
clad_th = 0.050
chan_gap = 0.075
r_add = 0.150
r_add_pos = 0.910

core_h = 42.40
refl_ax_th = 9.30
total_h = core_h + 2.0 * refl_ax_th # 61.00 cm
core_r = 19.45
refl_rad_th = 9.80
refl_rad_outer_r = core_r + refl_rad_th # 29.25 cm

ax_duct_r = 0.850 # Enhanced streaming duct
drum_count = 10
drum_r = 3.95
abs_th = 1.50
abs_angle = 148.0 # Enhanced drum authority
drum_clearance = 0.18
drum_rad_pos = core_r + drum_clearance + drum_r # 23.58 cm
b10_enrichment = 0.90

# Densities
rho_un = 13.50
rho_clad = 8.57
rho_yh2 = 4.30
rho_beo = 2.90
rho_b4c = 2.52
rho_he = 0.001146 # 2.5 MPa, 1050 K

# Zones: inner 19 pins (rings 0-2) at 19.20 wt%, outer 198 pins (rings 3-8) at 19.75 wt%
enr_inner = 19.20
enr_outer = 19.75


def build_materials(is_void=False, fuel_temp=1050.0, mod_temp=800.0, cool_temp=1050.0, mod_density=None):
    openmc.reset_auto_ids()
    
    # 1a. Inner UN Fuel
    f_in = openmc.Material(name="UN Fuel Inner (19.20%)")
    f_in.set_density("g/cm3", rho_un)
    u235_at_in = (enr_inner / 235.0) / (enr_inner / 235.0 + (100.0 - enr_inner) / 238.0)
    f_in.add_nuclide("U235", 0.5 * u235_at_in, "ao")
    f_in.add_nuclide("U238", 0.5 * (1.0 - u235_at_in), "ao")
    f_in.add_nuclide("N14", 0.5 * 0.9963, "ao")
    f_in.add_nuclide("N15", 0.5 * 0.0037, "ao")
    f_in.temperature = fuel_temp

    # 1b. Outer UN Fuel
    f_out = openmc.Material(name="UN Fuel Outer (19.75%)")
    f_out.set_density("g/cm3", rho_un)
    u235_at_out = (enr_outer / 235.0) / (enr_outer / 235.0 + (100.0 - enr_outer) / 238.0)
    f_out.add_nuclide("U235", 0.5 * u235_at_out, "ao")
    f_out.add_nuclide("U238", 0.5 * (1.0 - u235_at_out), "ao")
    f_out.add_nuclide("N14", 0.5 * 0.9963, "ao")
    f_out.add_nuclide("N15", 0.5 * 0.0037, "ao")
    f_out.temperature = fuel_temp

    # 2. Cladding
    clad = openmc.Material(name="Nb-1Zr Clad")
    clad.set_density("g/cm3", rho_clad)
    clad.add_element("Nb", 0.99, "wo")
    clad.add_element("Zr", 0.01, "wo")
    clad.temperature = fuel_temp

    # 3. Coolant
    cool = openmc.Material(name="Helium Coolant")
    c_rho = 1.0e-7 if is_void else rho_he
    cool.set_density("g/cm3", c_rho)
    cool.add_nuclide("He4", 1.0, "ao")
    cool.temperature = cool_temp

    # 4. Moderator
    m_rho = rho_yh2 if mod_density is None else mod_density
    mod = openmc.Material(name="Canned YH2.0")
    mod.set_density("g/cm3", m_rho)
    mod.add_nuclide("Y89", 1.0, "ao")
    mod.add_nuclide("H1", 2.0, "ao")
    mod.add_s_alpha_beta("c_H_in_ZrH")
    mod.temperature = mod_temp

    # 5. Reflector
    refl = openmc.Material(name="BeO Reflector")
    refl.set_density("g/cm3", rho_beo)
    refl.add_nuclide("Be9", 1.0, "ao")
    refl.add_nuclide("O16", 1.0, "ao")
    refl.add_s_alpha_beta("c_Be_in_BeO")
    refl.add_s_alpha_beta("c_O_in_BeO")
    refl.temperature = 800.0

    # 6. Absorber
    absorb = openmc.Material(name="B4C Absorber")
    absorb.set_density("g/cm3", rho_b4c)
    absorb.add_nuclide("B10", 4.0 * b10_enrichment, "ao")
    absorb.add_nuclide("B11", 4.0 * (1.0 - b10_enrichment), "ao")
    absorb.add_nuclide("C0", 1.0, "ao")
    absorb.temperature = 800.0

    return openmc.Materials([f_in, f_out, clad, cool, mod, refl, absorb]), f_in, f_out, clad, cool, mod, refl, absorb


def make_pin(f_mat, clad, cool, mod, refl):
    rc_out = rf + clad_th
    rch_out = rc_out + chan_gap
    z_cb = openmc.ZPlane(z0=-core_h/2.0)
    z_ct = openmc.ZPlane(z0=+core_h/2.0)
    z_rb = openmc.ZPlane(z0=-total_h/2.0, boundary_type="vacuum")
    z_rt = openmc.ZPlane(z0=+total_h/2.0, boundary_type="vacuum")

    cyl_f = openmc.ZCylinder(r=rf)
    cyl_c = openmc.ZCylinder(r=rc_out)
    cyl_ch = openmc.ZCylinder(r=rch_out)
    
    aux_cyls = []
    for k in range(6):
        ang = math.radians(k * 60.0)
        aux_cyls.append(openmc.ZCylinder(x0=r_add_pos*math.cos(ang), y0=r_add_pos*math.sin(ang), r=r_add))
    
    c_f = openmc.Cell(fill=f_mat, region=-cyl_f & +z_cb & -z_ct)
    c_c = openmc.Cell(fill=clad, region=+cyl_f & -cyl_c & +z_cb & -z_ct)
    c_ann = openmc.Cell(fill=cool, region=+cyl_c & -cyl_ch & +z_cb & -z_ct)
    
    aux_cells = []
    aux_union = None
    for c in aux_cyls:
        aux_cells.append(openmc.Cell(fill=cool, region=-c & +z_cb & -z_ct))
        aux_union = -c if aux_union is None else (aux_union | -c)
    
    mod_reg = +cyl_ch & +z_cb & -z_ct
    if aux_union is not None:
        mod_reg = mod_reg & ~aux_union
    c_m = openmc.Cell(fill=mod, region=mod_reg)
    
    cyl_ax = openmc.ZCylinder(r=ax_duct_r)
    c_db = openmc.Cell(fill=cool, region=-cyl_ax & +z_rb & -z_cb)
    c_dt = openmc.Cell(fill=cool, region=-cyl_ax & +z_ct & -z_rt)
    
    aux_db = []
    aux_dt = []
    aux_ax_union = None
    for c in aux_cyls:
        aux_db.append(openmc.Cell(fill=cool, region=-c & +z_rb & -z_cb))
        aux_dt.append(openmc.Cell(fill=cool, region=-c & +z_ct & -z_rt))
        aux_ax_union = -c if aux_ax_union is None else (aux_ax_union | -c)
        
    c_bb = openmc.Cell(fill=refl, region=+cyl_ax & +z_rb & -z_cb & ~aux_ax_union)
    c_bt = openmc.Cell(fill=refl, region=+cyl_ax & +z_ct & -z_rt & ~aux_ax_union)
    return openmc.Universe(cells=[c_f, c_c, c_ann, c_m, c_db, c_bb, c_dt, c_bt] + aux_cells + aux_db + aux_dt)


def build_geometry(f_in, f_out, clad, cool, mod, refl, absorb, drum_angles=None):
    if drum_angles is None:
        drum_angles = [0.0] * drum_count

    u_pin_in = make_pin(f_in, clad, cool, mod, refl)
    u_pin_out = make_pin(f_out, clad, cool, mod, refl)

    z_cb = openmc.ZPlane(z0=-core_h/2.0)
    z_ct = openmc.ZPlane(z0=+core_h/2.0)
    z_rb = openmc.ZPlane(z0=-total_h/2.0, boundary_type="vacuum")
    z_rt = openmc.ZPlane(z0=+total_h/2.0, boundary_type="vacuum")

    u_mod_buf = openmc.Universe(cells=[
        openmc.Cell(fill=mod, region=+z_cb & -z_ct),
        openmc.Cell(fill=refl, region=+z_rb & -z_cb),
        openmc.Cell(fill=refl, region=+z_ct & -z_rt)
    ])

    hex_lat = openmc.HexLattice()
    hex_lat.center = (0.0, 0.0)
    hex_lat.pitch = [pitch]
    hex_lat.orientation = "y"
    rings = []
    pin_centers = []
    v1 = np.array([pitch, 0.0])
    v2 = np.array([pitch / 2.0, pitch * np.sqrt(3.0) / 2.0])

    for r in range(8, -1, -1):
        u_ring = u_pin_in if r <= 2 else u_pin_out
        if r == 0:
            rings.append([u_ring])
            pin_centers.append((0.0, 0.0))
        else:
            rings.append([u_ring] * (6 * r))
            corners = [r*v2, r*(v2-v1), -r*v1, -r*v2, -r*(v2-v1), r*v1]
            for corn in range(6):
                p_start = corners[corn]
                p_end = corners[(corn + 1) % 6]
                for step in range(r):
                    pt = p_start + (p_end - p_start) * (step / r)
                    pin_centers.append((float(pt[0]), float(pt[1])))

    hex_lat.universes = rings
    hex_lat.outer = u_mod_buf

    cyl_core = openmc.ZCylinder(r=core_r)
    cyl_refl_outer = openmc.ZCylinder(r=refl_rad_outer_r, boundary_type="vacuum")

    cyl_drum_outer = openmc.ZCylinder(r=drum_r)
    cyl_drum_inner = openmc.ZCylinder(r=drum_r - abs_th)
    half_sec = np.radians(abs_angle / 2.0)
    plane_sec1 = openmc.Plane(a=np.sin(half_sec), b=-np.cos(half_sec), c=0.0, d=0.0)
    plane_sec2 = openmc.Plane(a=np.sin(half_sec), b=+np.cos(half_sec), c=0.0, d=0.0)
    plane_front = openmc.XPlane(x0=0.0)
    abs_reg = -cyl_drum_outer & +cyl_drum_inner & +plane_sec1 & +plane_sec2 & +plane_front
    beo_drum_reg = -cyl_drum_outer & ~abs_reg
    u_drum = openmc.Universe(cells=[openmc.Cell(fill=absorb, region=abs_reg), openmc.Cell(fill=refl, region=beo_drum_reg)])

    drum_cells = []
    drum_regs = []
    for k in range(drum_count):
        phi_k = k * 2.0 * np.pi / drum_count
        xk = drum_rad_pos * math.cos(phi_k)
        yk = drum_rad_pos * math.sin(phi_k)
        cyl_k = openmc.ZCylinder(x0=xk, y0=yk, r=drum_r)
        rot_deg = k * (360.0 / drum_count) + drum_angles[k]
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
    return openmc.Geometry(root), pin_centers


def compute_inventory(pin_centers):
    n_pins = len(pin_centers) # 217
    n_in = 19
    n_out = 198
    
    # Fuel & U-235
    v_f_pin = math.pi * rf**2 * core_h
    m_f_tot = n_pins * v_f_pin * rho_un * 1e-3
    m_u_in = n_in * v_f_pin * rho_un * 1e-3 * (238.0 / 252.0)
    m_u_out = n_out * v_f_pin * rho_un * 1e-3 * (238.0 / 252.0)
    m_u235 = m_u_in * (enr_inner / 100.0) + m_u_out * (enr_outer / 100.0)

    # Clad
    v_c_pin = math.pi * ((rf + clad_th)**2 - rf**2) * core_h
    m_clad = n_pins * v_c_pin * rho_clad * 1e-3

    # Coolant active & axial
    v_ann = math.pi * ((rf + clad_th + chan_gap)**2 - (rf + clad_th)**2) * core_h
    v_aux_core = 6 * math.pi * r_add**2 * core_h
    v_cool_core = n_pins * (v_ann + v_aux_core)

    v_duct_ax = math.pi * ax_duct_r**2 * (2.0 * refl_ax_th)
    v_aux_ax = 6 * math.pi * r_add**2 * (2.0 * refl_ax_th)
    v_cool_ax = n_pins * (v_duct_ax + v_aux_ax)

    # Moderator
    cell_area = (np.sqrt(3.0) / 2.0) * pitch**2
    v_cell_core = n_pins * cell_area * core_h
    v_mod = v_cell_core - (n_pins * v_f_pin + n_pins * v_c_pin + v_cool_core)
    m_mod = v_mod * rho_yh2 * 1e-3

    # Reflector BeO
    v_ax_gross = math.pi * refl_rad_outer_r**2 * (2.0 * refl_ax_th)
    v_beo_ax = v_ax_gross - v_cool_ax
    m_beo_ax = v_beo_ax * rho_beo * 1e-3

    v_cyl_active = math.pi * refl_rad_outer_r**2 * core_h
    v_refl_rad_gross = v_cyl_active - (math.pi * core_r**2 * core_h)
    v_drums = drum_count * math.pi * drum_r**2 * core_h
    v_refl_rad_stat = v_refl_rad_gross - v_drums

    v_abs_single = (abs_angle / 360.0) * math.pi * (drum_r**2 - (drum_r - abs_th)**2) * core_h
    v_abs_tot = drum_count * v_abs_single
    m_abs = v_abs_tot * rho_b4c * 1e-3

    v_beo_drums = v_drums - v_abs_tot
    m_beo_rad = (v_refl_rad_stat + v_beo_drums) * rho_beo * 1e-3
    m_beo_tot = m_beo_ax + m_beo_rad

    m_tot = m_f_tot + m_clad + m_mod + m_beo_tot + m_abs
    v_tot_L = math.pi * refl_rad_outer_r**2 * total_h * 1e-3

    return {
        "m_total_kg": m_tot,
        "v_total_L": v_tot_L,
        "m_u235_kg": m_u235,
        "m_fuel_kg": m_f_tot,
        "m_clad_kg": m_clad,
        "m_mod_kg": m_mod,
        "m_beo_tot_kg": m_beo_tot,
        "m_abs_kg": m_abs
    }


def main():
    print("=========================================================================")
    print("      EVALUATING CANDIDATE OPT-17 CHAMPION CANDIDATE FOR NO-REGRESSION   ")
    print("=========================================================================")
    
    mats_op, f_in, f_out, clad, cool, mod, refl, absorb = build_materials()
    geom_op, pin_centers = build_geometry(f_in, f_out, clad, cool, mod, refl, absorb, drum_angles=[0.0]*drum_count)
    inv = compute_inventory(pin_centers)

    print(f"\n[+] Total Mass: {inv['m_total_kg']:.2f} kg (Opt-16 Baseline: 540.44 kg | Margin: {540.44 - inv['m_total_kg']:+.2f} kg)")
    print(f"[+] Total Volume: {inv['v_total_L']:.2f} L (Opt-16: 163.96 L)")
    print(f"[+] Fissile U-235: {inv['m_u235_kg']:.2f} kg")

    def get_settings(particles=75000, batches=45, inactive=10):
        s = openmc.Settings()
        s.run_mode = "eigenvalue"
        s.particles = particles
        s.batches = batches
        s.inactive = inactive
        s.temperature = {"multipole": True, "method": "nearest", "tolerance": 1500.0}
        space = openmc.stats.Box((-core_r*0.6, -core_r*0.6, -core_h*0.3), (core_r*0.6, core_r*0.6, core_h*0.3))
        s.source = openmc.IndependentSource(space=space)
        return s

    # 1. Operating State & Power Peaking
    print("\n[1/4] Running Operating State Criticality & Peaking Factor...")
    mesh = openmc.RegularMesh()
    mesh.dimension = [80, 80]
    mesh.lower_left = [-refl_rad_outer_r, -refl_rad_outer_r]
    mesh.upper_right = [+refl_rad_outer_r, +refl_rad_outer_r]
    mesh_filter = openmc.MeshFilter(mesh)
    t_fis = openmc.Tally(name="fission_mesh")
    t_fis.filters = [mesh_filter]
    t_fis.scores = ["fission"]
    tallies = openmc.Tallies([t_fis])

    s_op = get_settings(particles=80000, batches=50, inactive=10)
    mod_op = openmc.Model(materials=mats_op, geometry=geom_op, settings=s_op, tallies=tallies)
    mod_op.export_to_xml()
    openmc.run(openmc_exec=OPENMC_EXEC, output=False)
    with openmc.StatePoint(f"statepoint.{s_op.batches}.h5") as sp:
        k_op = float(sp.keff.nominal_value)
        std_op = float(sp.keff.std_dev)
        tf = sp.get_tally(name="fission_mesh")
        fiss_data = tf.get_values(scores=["fission"]).reshape((80, 80))
    os.remove(f"statepoint.{s_op.batches}.h5")

    xs = np.linspace(-refl_rad_outer_r, refl_rad_outer_r, 80)
    ys = np.linspace(-refl_rad_outer_r, refl_rad_outer_r, 80)
    XX, YY = np.meshgrid(xs, ys)
    core_mask = (XX**2 + YY**2) <= (core_r**2)
    fiss_core = fiss_data[core_mask]
    avg_fiss = np.mean(fiss_core[fiss_core > 0])
    max_fiss = np.max(fiss_core)
    Fq = max_fiss / avg_fiss
    excess_pcm = (k_op - 1.0) / k_op * 1e5
    life_yr = excess_pcm / 68.2

    print(f"  --> k_op = {k_op:.5f} +/- {std_op:.5f}")
    print(f"  --> Excess Reactivity = {excess_pcm:+.1f} pcm")
    print(f"  --> Depletion Life = {life_yr:.1f} YEARS (Opt-16 Baseline: 28.0 yr | PASS: {life_yr >= 28.0})")
    print(f"  --> Radial Peaking Factor Fq = {Fq:.3f} (Opt-16 Baseline: 3.428 | IMPROVEMENT: {3.428 - Fq:+.3f})")

    # 2. Shutdown State (180° drums)
    print("\n[2/4] Running All-Drums-In Shutdown State...")
    geom_in, _ = build_geometry(f_in, f_out, clad, cool, mod, refl, absorb, drum_angles=[180.0]*drum_count)
    s_in = get_settings(particles=60000, batches=40, inactive=10)
    mod_in = openmc.Model(materials=mats_op, geometry=geom_in, settings=s_in)
    mod_in.export_to_xml()
    openmc.run(openmc_exec=OPENMC_EXEC, output=False)
    with openmc.StatePoint(f"statepoint.{s_in.batches}.h5") as sp:
        k_in = float(sp.keff.nominal_value)
        std_in = float(sp.keff.std_dev)
    os.remove(f"statepoint.{s_in.batches}.h5")

    worth_pcm = (1.0/k_in - 1.0/k_op) * 1e5
    sdm_pcm = (1.0 - k_in) / k_in * 1e5
    print(f"  --> k_in = {k_in:.5f} +/- {std_in:.5f}")
    print(f"  --> Drum Worth = {worth_pcm:.1f} pcm (Opt-16 Baseline: 4740.9 pcm | Margin: {worth_pcm - 4740.9:+.1f} pcm)")
    print(f"  --> All-In SDM = {sdm_pcm:.1f} pcm (Opt-16 Baseline: 2831.8 pcm | Margin: {sdm_pcm - 2831.8:+.1f} pcm)")

    # 3. Single Stuck-Drum
    print("\n[3/4] Running Single Stuck-Drum State (Drum 0 OUT, Drums 1-9 IN)...")
    stk_angles = [0.0] + [180.0] * (drum_count - 1)
    geom_stk, _ = build_geometry(f_in, f_out, clad, cool, mod, refl, absorb, drum_angles=stk_angles)
    mod_stk = openmc.Model(materials=mats_op, geometry=geom_stk, settings=s_in)
    mod_stk.export_to_xml()
    openmc.run(openmc_exec=OPENMC_EXEC, output=False)
    with openmc.StatePoint(f"statepoint.{s_in.batches}.h5") as sp:
        k_stk = float(sp.keff.nominal_value)
        std_stk = float(sp.keff.std_dev)
    os.remove(f"statepoint.{s_in.batches}.h5")

    sdm_stk_pcm = (1.0 - k_stk) / k_stk * 1e5
    print(f"  --> k_stk = {k_stk:.5f} +/- {std_stk:.5f}")
    print(f"  --> Stuck SDM = {sdm_stk_pcm:.1f} pcm (Opt-16 Baseline: 2393.7 pcm | Margin: {sdm_stk_pcm - 2393.7:+.1f} pcm)")

    # 4. Voided Depressurization State (High-precision)
    print("\n[4/4] Running Coolant Void Depressurization State (Vacuum)...")
    mats_vd, f_in_vd, f_out_vd, clad_vd, cool_vd, mod_vd, refl_vd, abs_vd = build_materials(is_void=True)
    geom_vd, _ = build_geometry(f_in_vd, f_out_vd, clad_vd, cool_vd, mod_vd, refl_vd, abs_vd, drum_angles=[0.0]*drum_count)
    s_vd = get_settings(particles=150000, batches=60, inactive=15)
    mod_vd = openmc.Model(materials=mats_vd, geometry=geom_vd, settings=s_vd)
    mod_vd.export_to_xml()
    openmc.run(openmc_exec=OPENMC_EXEC, output=False)
    with openmc.StatePoint(f"statepoint.{s_vd.batches}.h5") as sp:
        k_vd = float(sp.keff.nominal_value)
        std_vd = float(sp.keff.std_dev)
    os.remove(f"statepoint.{s_vd.batches}.h5")

    # High-precision operating reference for void
    mod_op_hi = openmc.Model(materials=mats_op, geometry=geom_op, settings=s_vd)
    mod_op_hi.export_to_xml()
    openmc.run(openmc_exec=OPENMC_EXEC, output=False)
    with openmc.StatePoint(f"statepoint.{s_vd.batches}.h5") as sp:
        k_op_hi = float(sp.keff.nominal_value)
        std_op_hi = float(sp.keff.std_dev)
    os.remove(f"statepoint.{s_vd.batches}.h5")

    drho_vd = (1.0/k_op_hi - 1.0/k_vd) * 1e5
    u_drho_vd = math.sqrt((std_op_hi/k_op_hi**2)**2 + (std_vd/k_vd**2)**2) * 1e5
    b1_vd = drho_vd + 1.0 * u_drho_vd
    b2_vd = drho_vd + 2.0 * u_drho_vd
    b3_vd = drho_vd + 3.0 * u_drho_vd

    print(f"  Operating (2.5 MPa He): k_op = {k_op_hi:.5f} +/- {std_op_hi:.5f}")
    print(f"  Voided (0.0 MPa Vac):   k_vd = {k_vd:.5f} +/- {std_vd:.5f}")
    print(f"  --> Delta_rho_void = {drho_vd:+.1f} +/- {u_drho_vd:.1f} pcm")
    print(f"      1-sigma Upper Bound: {b1_vd:+.1f} pcm")
    print(f"      2-sigma Upper Bound: {b2_vd:+.1f} pcm")
    print(f"      3-sigma Upper Bound: {b3_vd:+.1f} pcm")
    print(f"      (Opt-16 Baseline: -87.1 +/- 37.2 pcm)")

    print("\n=========================================================================")
    print("                           SUMMARY OF CHECKPOINTS                        ")
    print("=========================================================================")
    print(f"1. Total Mass:       {inv['m_total_kg']:.2f} kg <= 540.44 kg  --> {'PASS' if inv['m_total_kg'] <= 540.44 else 'FAIL'}")
    print(f"2. Lifetime:         {life_yr:.1f} yr >= 28.0 yr             --> {'PASS' if life_yr >= 28.0 else 'FAIL'}")
    print(f"3. Drum Worth:       {worth_pcm:.1f} pcm >= 4740.9 pcm       --> {'PASS' if worth_pcm >= 4740.9 else 'FAIL'}")
    print(f"4. Shutdown Margin:  {sdm_pcm:.1f} pcm >= 2831.8 pcm         --> {'PASS' if sdm_pcm >= 2831.8 else 'FAIL'}")
    print(f"5. Stuck-Drum SDM:   {sdm_stk_pcm:.1f} pcm >= 2393.7 pcm     --> {'PASS' if sdm_stk_pcm >= 2393.7 else 'FAIL'}")
    print(f"6. Peaking Factor:   Fq = {Fq:.3f} <= 3.428                 --> {'PASS' if Fq <= 3.428 else 'FAIL'}")
    print(f"7. Void Reactivity:  {drho_vd:+.1f} pcm <= -87.1 pcm        --> {'PASS' if drho_vd <= -87.1 else 'FAIL'}")


if __name__ == "__main__":
    main()
