#!/usr/bin/env python3
"""
screen_neutronics.py
Neutronic screening of candidate geometries using OpenMC.
Evaluates keff(0 deg) and keff(180 deg) at BOL.
"""

import os
import sys
import time
import math
import numpy as np
import openmc

OPENMC_EXEC = "/home/aleri/miniforge3/envs/openmc-env/bin/openmc"
CROSS_SECTIONS_PATH = "/home/aleri/wmp_data/cross_sections_wmp.xml"

# Densities
RHO_UN = 13.50
RHO_CLAD = 8.57
RHO_YH2 = 4.30
RHO_BEO = 2.90
RHO_B4C = 2.52
RHO_AL2O3 = 3.95
P_NOM = 2.50e6
R_GAS = 8314.462618
M_HE = 4.002602
N_PINS = 217
DRUM_COUNT = 10

def rho_helium(P, T):
    return (P * M_HE) / (R_GAS * max(T, 100.0))

def build_flexible_model(run_dir,
                         r_fuel=0.580,
                         clad_th=0.048,
                         chan_gap=0.084,
                         al2o3_th=0.010,
                         pitch=2.270,
                         H_core=42.45,
                         refl_ax_th=9.275,
                         refl_rad_th=9.80,
                         enr_inner=19.30,
                         enr_outer=19.75,
                         drum_r=3.95,
                         abs_th=1.50,
                         abs_angle=148.0,
                         drum_clearance=0.180,
                         drum_angle=0.0,
                         particles=12000,
                         batches=25,
                         inactive=8,
                         seed=42):
    os.makedirs(run_dir, exist_ok=True)
    openmc.reset_auto_ids()
    openmc.config['cross_sections'] = CROSS_SECTIONS_PATH

    core_h = H_core
    drum_angles = [drum_angle] * DRUM_COUNT
    r_f = r_fuel
    r_c = r_f + clad_th
    r_ch = r_c + chan_gap
    r_m_in = r_ch + al2o3_th

    r_add = 0.150
    r_add_pos = 0.910 * (pitch / 2.270)
    core_r = 19.45 * (pitch / 2.270)
    refl_rad_outer_r = core_r + refl_rad_th
    total_h = core_h + 2.0 * refl_ax_th
    ax_duct_r = 0.860
    drum_rad_pos = core_r + drum_clearance + drum_r

    T_fuel = 1215.0
    T_clad = 1190.0
    T_cool = 1025.0
    T_mod = 940.0
    T_refl = 880.0

    # Materials
    f_in = openmc.Material(name="UN Fuel Inner")
    f_in.set_density("g/cm3", RHO_UN)
    u235_at_in = (enr_inner / 235.0) / (enr_inner / 235.0 + (100.0 - enr_inner) / 238.0)
    f_in.add_nuclide("U235", 0.5 * u235_at_in, "ao")
    f_in.add_nuclide("U238", 0.5 * (1.0 - u235_at_in), "ao")
    f_in.add_nuclide("N14", 0.5 * 0.9963, "ao")
    f_in.add_nuclide("N15", 0.5 * 0.0037, "ao")
    f_in.temperature = T_fuel

    f_out = openmc.Material(name="UN Fuel Outer")
    f_out.set_density("g/cm3", RHO_UN)
    u235_at_out = (enr_outer / 235.0) / (enr_outer / 235.0 + (100.0 - enr_outer) / 238.0)
    f_out.add_nuclide("U235", 0.5 * u235_at_out, "ao")
    f_out.add_nuclide("U238", 0.5 * (1.0 - u235_at_out), "ao")
    f_out.add_nuclide("N14", 0.5 * 0.9963, "ao")
    f_out.add_nuclide("N15", 0.5 * 0.0037, "ao")
    f_out.temperature = T_fuel

    clad = openmc.Material(name="Nb-1Zr Clad")
    clad.set_density("g/cm3", RHO_CLAD)
    clad.add_element("Nb", 0.99, "wo")
    clad.add_element("Zr", 0.01, "wo")
    clad.temperature = T_clad

    al2o3 = openmc.Material(name="Al2O3 Insulator Barrier")
    al2o3.set_density("g/cm3", RHO_AL2O3)
    al2o3.add_element("Al", 2.0, "ao")
    al2o3.add_element("O", 3.0, "ao")
    al2o3.temperature = T_mod

    rho_he_si = rho_helium(P_NOM, T_cool) / 1000.0
    cool_mat = openmc.Material(name="He Coolant")
    cool_mat.set_density("g/cm3", max(rho_he_si, 1e-7))
    cool_mat.add_nuclide("He4", 1.0, "ao")
    cool_mat.temperature = T_cool

    mod = openmc.Material(name="Canned YH2.0")
    mod.set_density("g/cm3", RHO_YH2)
    mod.add_nuclide("Y89", 1.0, "ao")
    mod.add_nuclide("H1", 2.0, "ao")
    mod.add_s_alpha_beta("c_H_in_ZrH")
    mod.temperature = T_mod

    refl = openmc.Material(name="BeO Reflector")
    refl.set_density("g/cm3", RHO_BEO)
    refl.add_nuclide("Be9", 1.0, "ao")
    refl.add_nuclide("O16", 1.0, "ao")
    refl.add_s_alpha_beta("c_Be_in_BeO")
    refl.add_s_alpha_beta("c_O_in_BeO")
    refl.temperature = T_refl

    absorb = openmc.Material(name="B4C Absorber")
    absorb.set_density("g/cm3", RHO_B4C)
    absorb.add_nuclide("B10", 4.0 * 0.90, "ao")
    absorb.add_nuclide("B11", 4.0 * 0.10, "ao")
    absorb.add_nuclide("C0", 1.0, "ao")
    absorb.temperature = T_refl

    materials = openmc.Materials([f_in, f_out, clad, al2o3, cool_mat, mod, refl, absorb])

    # Geometry
    z_cb = openmc.ZPlane(z0=-core_h / 2.0)
    z_ct = openmc.ZPlane(z0=+core_h / 2.0)
    z_rb = openmc.ZPlane(z0=-total_h / 2.0, boundary_type="vacuum")
    z_rt = openmc.ZPlane(z0=+total_h / 2.0, boundary_type="vacuum")

    cyl_f = openmc.ZCylinder(r=r_f)
    cyl_c = openmc.ZCylinder(r=r_c)
    cyl_ch = openmc.ZCylinder(r=r_ch)
    cyl_m_in = openmc.ZCylinder(r=r_m_in)

    aux_cyls = []
    for k in range(6):
        ang = math.radians(k * 60.0)
        aux_cyls.append(openmc.ZCylinder(x0=r_add_pos * math.cos(ang),
                                         y0=r_add_pos * math.sin(ang),
                                         r=r_add))
    cyl_ax = openmc.ZCylinder(r=ax_duct_r)

    def make_pin(f_mat):
        c_f = openmc.Cell(fill=f_mat, region=-cyl_f & +z_cb & -z_ct)
        c_f.temperature = T_fuel
        c_c = openmc.Cell(fill=clad, region=+cyl_f & -cyl_c & +z_cb & -z_ct)
        c_c.temperature = T_clad
        c_ann = openmc.Cell(fill=cool_mat, region=+cyl_c & -cyl_ch & +z_cb & -z_ct)
        c_ann.temperature = T_cool
        c_al = openmc.Cell(fill=al2o3, region=+cyl_ch & -cyl_m_in & +z_cb & -z_ct)
        c_al.temperature = T_mod

        aux_cells = []
        aux_union = None
        for c in aux_cyls:
            aux_cells.append(openmc.Cell(fill=cool_mat, region=-c & +z_cb & -z_ct))
            aux_union = -c if aux_union is None else (aux_union | -c)

        mod_reg = +cyl_m_in & +z_cb & -z_ct
        if aux_union is not None:
            mod_reg = mod_reg & ~aux_union
        c_m = openmc.Cell(fill=mod, region=mod_reg)

        c_db = openmc.Cell(fill=cool_mat, region=-cyl_ax & +z_rb & -z_cb)
        c_dt = openmc.Cell(fill=cool_mat, region=-cyl_ax & +z_ct & -z_rt)

        aux_db = []
        aux_dt = []
        aux_ax_union = None
        for c in aux_cyls:
            aux_db.append(openmc.Cell(fill=cool_mat, region=-c & +z_rb & -z_cb))
            aux_dt.append(openmc.Cell(fill=cool_mat, region=-c & +z_ct & -z_rt))
            aux_ax_union = -c if aux_ax_union is None else (aux_ax_union | -c)

        c_bb = openmc.Cell(fill=refl, region=+cyl_ax & +z_rb & -z_cb & ~aux_ax_union)
        c_bb.temperature = T_refl
        c_bt = openmc.Cell(fill=refl, region=+cyl_ax & +z_ct & -z_rt & ~aux_ax_union)
        c_bt.temperature = T_refl

        return openmc.Universe(cells=[c_f, c_c, c_ann, c_al, c_m, c_db, c_bb, c_dt, c_bt] + aux_cells + aux_db + aux_dt)

    u_in = make_pin(f_in)
    u_out = make_pin(f_out)

    v1 = np.array([pitch, 0.0])
    v2 = np.array([pitch / 2.0, pitch * np.sqrt(3.0) / 2.0])

    rings = []
    for r in range(8, -1, -1):
        u_p = u_in if r <= 2 else u_out
        ring_u = []
        if r == 0:
            ring_u.append(u_p)
        else:
            for corn in range(6):
                for step in range(r):
                    ring_u.append(u_p)
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

    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    settings.temperature = {"multipole": True, "method": "nearest", "tolerance": 1500.0}
    settings.seed = seed

    space = openmc.stats.Box((-core_r * 0.6, -core_r * 0.6, -core_h * 0.3),
                            (core_r * 0.6, core_r * 0.6, core_h * 0.3))
    settings.source = openmc.IndependentSource(space=space)

    materials.export_to_xml(os.path.join(run_dir, "materials.xml"))
    geometry.export_to_xml(os.path.join(run_dir, "geometry.xml"))
    settings.export_to_xml(os.path.join(run_dir, "settings.xml"))

def run_openmc_calc(run_dir):
    os.environ["OPENMC_CROSS_SECTIONS"] = CROSS_SECTIONS_PATH
    os.environ["OMP_NUM_THREADS"] = "8"
    cwd = os.getcwd()
    os.chdir(run_dir)
    try:
        openmc.run(openmc_exec=OPENMC_EXEC, output=False)
    finally:
        os.chdir(cwd)
    
    # Read statepoint
    sp_files = [f for f in os.listdir(run_dir) if f.startswith("statepoint.") and f.endswith(".h5")]
    sp_path = os.path.join(run_dir, sorted(sp_files)[-1])
    with openmc.StatePoint(sp_path) as sp:
        k_val = float(sp.keff.nominal_value)
        k_sig = float(sp.keff.std_dev)
    return k_val, k_sig

if __name__ == "__main__":
    from screen_candidates import candidates, calc_geometry_properties

    print("=" * 80)
    print("SCREENING NEUTRONICS: 4 CANDIDATES + SIM 4B BASELINE")
    print("=" * 80)

    for name, params in candidates.items():
        print(f"\n--- Testing Candidate: {name} ---")
        geo_info = calc_geometry_properties(**params)
        print(f"    Mass: {geo_info['M_tot']:.2f} kg | U-235: {geo_info['u235_kg']:.2f} kg | DP: {geo_info['dp_kpa']:.2f} kPa")
        
        # 1. theta = 0 deg
        dir_0 = f"/tmp/screen_{name}_0deg"
        build_flexible_model(dir_0, drum_angle=0.0, **params)
        k_0, sig_0 = run_openmc_calc(dir_0)
        rho_ex = (k_0 - 1.0) / k_0 * 1e5
        
        # 2. theta = 180 deg
        dir_180 = f"/tmp/screen_{name}_180deg"
        build_flexible_model(dir_180, drum_angle=180.0, **params)
        k_180, sig_180 = run_openmc_calc(dir_180)
        sdm = (1.0 - k_180) / k_180 * 1e5
        worth = (1.0 / k_180 - 1.0 / k_0) * 1e5
        
        # Burnup rate: P_th = 400 kW, U-235 consumption ~ 0.1845 kg/yr
        # Specific burnup slope (pcm/yr) scales inversely with U-235 mass:
        burn_slope = 76.5 * (24.47 / geo_info['u235_kg'])
        est_life = rho_ex / burn_slope
        
        print(f"    k_eff(0 deg):   {k_0:.5f} +/- {sig_0:.5f} (Excess rho = {rho_ex:+.1f} pcm)")
        print(f"    k_eff(180 deg): {k_180:.5f} +/- {sig_180:.5f} (SDM = {sdm:.1f} pcm, Worth = {worth:.1f} pcm)")
        print(f"    Burnup slope:   {burn_slope:.2f} pcm/yr")
        print(f"    Estimated Lifespan: {est_life:.2f} years")

