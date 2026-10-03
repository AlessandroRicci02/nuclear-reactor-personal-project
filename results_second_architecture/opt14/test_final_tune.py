#!/usr/bin/env python3
"""
Test Candidate Opt-14 with P = 2.27 cm, r_add = 0.160 cm, r_duct = 0.825 cm
to confirm k_op ~ 1.018 - 1.021 and Delta_rho_void ~ -90 pcm.
"""

import os
import sys
import math
import numpy as np
import openmc

os.environ["PATH"] = "/home/aleri/miniforge3/envs/openmc-env/bin:" + os.environ.get("PATH", "")
os.environ["OPENMC_CROSS_SECTIONS"] = "/home/aleri/wmp_data/cross_sections_wmp.xml"
os.environ["OMP_NUM_THREADS"] = "8"
openmc.config['cross_sections'] = "/home/aleri/wmp_data/cross_sections_wmp.xml"

workdir = "/home/aleri/OpenMC_projects/first_test/opt14/test_opt14_final_tune"
os.makedirs(workdir, exist_ok=True)
os.chdir(workdir)

def test_final_tune(pitch=2.27, r_add=0.160, r_duct=0.825, particles=80000, batches=45):
    openmc.reset_auto_ids()
    print(f"\n--- Testing Final Tune: P={pitch} cm, r_add={r_add} cm, r_duct={r_duct} cm ---")

    fuel = openmc.Material(name="UN Fuel")
    fuel.set_density("g/cm3", 13.50)
    fuel.add_nuclide("U235", 0.5 * 0.199522, "ao")
    fuel.add_nuclide("U238", 0.5 * 0.800478, "ao")
    fuel.add_nuclide("N14", 0.5 * 0.9963, "ao")
    fuel.add_nuclide("N15", 0.5 * 0.0037, "ao")
    fuel.temperature = 1050.0

    clad = openmc.Material(name="Nb-1Zr Clad")
    clad.set_density("g/cm3", 8.57)
    clad.add_element("Nb", 0.99, "wo")
    clad.add_element("Zr", 0.01, "wo")
    clad.temperature = 1050.0

    cool_op = openmc.Material(name="Pure He Coolant")
    cool_op.set_density("g/cm3", 0.001146)
    cool_op.add_nuclide("He4", 1.0, "ao")
    cool_op.temperature = 1050.0

    cool_vd = openmc.Material(name="Void Coolant")
    cool_vd.set_density("g/cm3", 1e-7)
    cool_vd.add_nuclide("He4", 1.0, "ao")
    cool_vd.temperature = 1050.0

    mod = openmc.Material(name="Canned YH2.0")
    mod.set_density("g/cm3", 4.30)
    mod.add_nuclide("Y89", 1.0, "ao")
    mod.add_nuclide("H1", 2.00, "ao")
    mod.add_s_alpha_beta("c_H_in_ZrH")
    mod.temperature = 800.0

    refl = openmc.Material(name="BeO Reflector")
    refl.set_density("g/cm3", 2.90)
    refl.add_nuclide("Be9", 1.0, "ao")
    refl.add_nuclide("O16", 1.0, "ao")
    refl.add_s_alpha_beta("c_Be_in_BeO")
    refl.add_s_alpha_beta("c_O_in_BeO")
    refl.temperature = 800.0

    abs_mat = openmc.Material(name="B4C Absorber")
    abs_mat.set_density("g/cm3", 2.52)
    abs_mat.add_nuclide("B10", 4.0 * 0.90, "ao")
    abs_mat.add_nuclide("B11", 4.0 * 0.10, "ao")
    abs_mat.add_nuclide("C0", 1.0, "ao")
    abs_mat.temperature = 800.0

    H_core = 42.0
    R_core = 19.45
    t_refl_rad = 9.80
    t_refl_ax = 9.50
    R_outer = R_core + t_refl_rad
    H_tot = H_core + 2.0 * t_refl_ax

    r_fuel = 0.580
    r_clad = 0.630
    r_chan = 0.705
    r_add_pos = 0.910

    def build_geom(cool, drum_angles=None):
        if drum_angles is None:
            drum_angles = [0.0] * 10

        s_f = openmc.ZCylinder(r=r_fuel)
        s_cl = openmc.ZCylinder(r=r_clad)
        s_ch = openmc.ZCylinder(r=r_chan)

        z_bot = openmc.ZPlane(z0=-H_core/2.0)
        z_top = openmc.ZPlane(z0=+H_core/2.0)
        z_refl_bot = openmc.ZPlane(z0=-H_tot/2.0, boundary_type="vacuum")
        z_refl_top = openmc.ZPlane(z0=+H_tot/2.0, boundary_type="vacuum")

        add_cyls = []
        for k in range(6):
            ang = math.radians(0.0 + k * 60.0)
            xv = r_add_pos * math.cos(ang)
            yv = r_add_pos * math.sin(ang)
            add_cyls.append(openmc.ZCylinder(x0=xv, y0=yv, r=r_add))

        c_f = openmc.Cell(fill=fuel, region=-s_f & +z_bot & -z_top)
        c_cl = openmc.Cell(fill=clad, region=+s_f & -s_cl & +z_bot & -z_top)
        c_co = openmc.Cell(fill=cool, region=+s_cl & -s_ch & +z_bot & -z_top)

        add_cells = []
        add_union = None
        for cyl in add_cyls:
            reg_c = -cyl & +z_bot & -z_top
            add_cells.append(openmc.Cell(fill=cool, region=reg_c))
            add_union = -cyl if add_union is None else (add_union | -cyl)

        mod_reg = +s_ch & +z_bot & -z_top
        if add_union is not None:
            mod_reg = mod_reg & ~add_union
        c_mo = openmc.Cell(fill=mod, region=mod_reg)

        # Axial reflector ducts
        s_duct = openmc.ZCylinder(r=r_duct)
        c_co_bot = openmc.Cell(fill=cool, region=-s_duct & +z_refl_bot & -z_bot)
        c_co_top = openmc.Cell(fill=cool, region=-s_duct & +z_top & -z_refl_top)
        
        # Penetrating auxiliary channels through axial reflectors
        aux_bot_cells = []
        aux_top_cells = []
        aux_ax_union = None
        for cyl in add_cyls:
            aux_bot_cells.append(openmc.Cell(fill=cool, region=-cyl & +z_refl_bot & -z_bot))
            aux_top_cells.append(openmc.Cell(fill=cool, region=-cyl & +z_top & -z_refl_top))
            aux_ax_union = -cyl if aux_ax_union is None else (aux_ax_union | -cyl)
        
        c_be_bot = openmc.Cell(fill=refl, region=+s_duct & +z_refl_bot & -z_bot & ~aux_ax_union)
        c_be_top = openmc.Cell(fill=refl, region=+s_duct & +z_top & -z_refl_top & ~aux_ax_union)
        u_pin = openmc.Universe(cells=[c_f, c_cl, c_co, c_mo, c_co_bot, c_be_bot, c_co_top, c_be_top] + add_cells + aux_bot_cells + aux_top_cells)

        c_mo_only = openmc.Cell(fill=mod, region=+z_bot & -z_top)
        c_be_bot_only = openmc.Cell(fill=refl, region=+z_refl_bot & -z_bot)
        c_be_top_only = openmc.Cell(fill=refl, region=+z_top & -z_refl_top)
        u_mod = openmc.Universe(cells=[c_mo_only, c_be_bot_only, c_be_top_only])

        n_rings = int(np.ceil(R_core / pitch)) + 2
        hex_lat = openmc.HexLattice()
        hex_lat.center = (0.0, 0.0)
        hex_lat.pitch = [pitch]
        hex_lat.orientation = 'y'

        rings = []
        pin_count = 0
        for r in range(n_rings - 1, -1, -1):
            if r == 0:
                rings.append([u_pin])
                pin_count += 1
            else:
                ring_u = []
                for k in range(6 * r):
                    side = k // r
                    pos_in_side = k % r
                    dist = r * pitch * np.sqrt(1.0 - 0.25 * (pos_in_side / r - 0.5)**2)
                    if dist <= (R_core - r_chan):
                        ring_u.append(u_pin)
                        pin_count += 1
                    else:
                        ring_u.append(u_mod)
                rings.append(ring_u)

        hex_lat.universes = rings
        hex_lat.outer = u_mod

        s_core = openmc.ZCylinder(r=R_core)
        s_outer = openmc.ZCylinder(r=R_outer, boundary_type="vacuum")

        cyl_drum_out = openmc.ZCylinder(r=3.95)
        cyl_drum_in = openmc.ZCylinder(r=3.95 - 1.45)
        half_phi = np.radians(140.0 / 2.0)
        p1 = openmc.Plane(a=np.sin(half_phi), b=-np.cos(half_phi), c=0.0, d=0.0)
        p2 = openmc.Plane(a=np.sin(half_phi), b=+np.cos(half_phi), c=0.0, d=0.0)
        px = openmc.XPlane(x0=0.0)

        abs_reg = -cyl_drum_out & +cyl_drum_in & +p1 & +p2 & +px
        beo_reg = -cyl_drum_out & ~abs_reg
        u_drum = openmc.Universe(cells=[openmc.Cell(fill=abs_mat, region=abs_reg), openmc.Cell(fill=refl, region=beo_reg)])

        drum_cells = []
        drum_regions = []
        r_drum_pos = R_core + 0.20 + 3.95
        for k in range(10):
            ang = k * 2.0 * np.pi / 10.0
            xk = r_drum_pos * np.cos(ang)
            yk = r_drum_pos * np.sin(ang)
            dk = openmc.ZCylinder(x0=xk, y0=yk, r=3.95)
            rot = k * 36.0 + drum_angles[k]
            cdk = openmc.Cell(fill=u_drum, region=-dk & +z_bot & -z_top)
            cdk.translation = (xk, yk, 0.0)
            cdk.rotation = (0.0, 0.0, rot)
            drum_cells.append(cdk)
            drum_regions.append(-dk)

        drum_union = drum_regions[0]
        for reg in drum_regions[1:]:
            drum_union = drum_union | reg

        c_core = openmc.Cell(fill=hex_lat, region=-s_core & +z_refl_bot & -z_refl_top)
        c_refl_rad = openmc.Cell(fill=refl, region=+s_core & -s_outer & +z_bot & -z_top & ~drum_union)
        c_refl_bot = openmc.Cell(fill=refl, region=+s_core & -s_outer & +z_refl_bot & -z_bot)
        c_refl_top = openmc.Cell(fill=refl, region=+s_core & -s_outer & +z_top & -z_refl_top)

        root_u = openmc.Universe(cells=[c_core, c_refl_rad, c_refl_bot, c_refl_top] + drum_cells)
        return openmc.Geometry(root_u), pin_count

    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = particles
    settings.batches = batches
    settings.inactive = 10
    settings.temperature = {'multipole': True, 'method': 'nearest', 'tolerance': 1500.0}
    space = openmc.stats.Box((-R_core*0.6, -R_core*0.6, -H_core*0.3), (R_core*0.6, R_core*0.6, H_core*0.3))
    settings.source = openmc.IndependentSource(space=space)

    # 1. Operating state (0 deg)
    mats_op = openmc.Materials([fuel, clad, cool_op, mod, refl, abs_mat])
    mats_op.cross_sections = "/home/aleri/wmp_data/cross_sections_wmp.xml"
    geom_op, pin_count = build_geom(cool_op, drum_angles=[0.0]*10)
    print(f"Total Pins: {pin_count}")
    mod_op = openmc.Model(materials=mats_op, geometry=geom_op, settings=settings)
    mod_op.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{batches}.h5") as sp:
        k_op = float(sp.keff.nominal_value)
        s_op = float(sp.keff.std_dev)

    # 2. Shutdown state (180 deg)
    geom_in, _ = build_geom(cool_op, drum_angles=[180.0]*10)
    mod_in = openmc.Model(materials=mats_op, geometry=geom_in, settings=settings)
    mod_in.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{batches}.h5") as sp:
        k_in = float(sp.keff.nominal_value)
        s_in = float(sp.keff.std_dev)

    # 3. Stuck drum (0 deg drum 0, 180 deg drums 1-9)
    geom_stk, _ = build_geom(cool_op, drum_angles=[0.0] + [180.0]*9)
    mod_stk = openmc.Model(materials=mats_op, geometry=geom_stk, settings=settings)
    mod_stk.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{batches}.h5") as sp:
        k_stk = float(sp.keff.nominal_value)
        s_stk = float(sp.keff.std_dev)

    # 4. Voided state (0.0 MPa vacuum)
    mats_vd = openmc.Materials([fuel, clad, cool_vd, mod, refl, abs_mat])
    mats_vd.cross_sections = "/home/aleri/wmp_data/cross_sections_wmp.xml"
    geom_vd, _ = build_geom(cool_vd, drum_angles=[0.0]*10)
    mod_vd = openmc.Model(materials=mats_vd, geometry=geom_vd, settings=settings)
    mod_vd.export_to_xml()
    openmc.run(output=False)
    with openmc.StatePoint(f"statepoint.{batches}.h5") as sp:
        k_vd = float(sp.keff.nominal_value)
        s_vd = float(sp.keff.std_dev)

    drum_worth = (1.0/k_in - 1.0/k_op) * 1e5
    sdm = (1.0 - k_in)/k_in * 1e5
    sdm_stk = (1.0 - k_stk)/k_stk * 1e5
    drho = (1.0/k_op - 1.0/k_vd) * 1e5
    u_drho = math.sqrt((s_op/k_op**2)**2 + (s_vd/k_vd**2)**2) * 1e5
    b1 = drho + u_drho
    b2 = drho + 2.0 * u_drho
    b3 = drho + 3.0 * u_drho

    print(f"Results for Final Tune (P={pitch} cm, r_add={r_add} cm):")
    print(f"  k_op       = {k_op:.5f} +/- {s_op:.5f} (Excess = {(k_op-1)/k_op*1e5:+.1f} pcm)")
    print(f"  k_in       = {k_in:.5f} +/- {s_in:.5f}")
    print(f"  k_stk      = {k_stk:.5f} +/- {s_stk:.5f}")
    print(f"  Drum Worth = {drum_worth:.1f} pcm")
    print(f"  SDM        = {sdm:.1f} pcm")
    print(f"  SDM_stuck  = {sdm_stk:.1f} pcm")
    print(f"  k_vd       = {k_vd:.5f} +/- {s_vd:.5f}")
    print(f"  Delta_rho_void = {drho:+.1f} +/- {u_drho:.1f} pcm")
    print(f"    1-sigma bound: {b1:+.1f} pcm")
    print(f"    2-sigma bound: {b2:+.1f} pcm")
    print(f"    3-sigma bound: {b3:+.1f} pcm")

if __name__ == "__main__":
    test_final_tune(pitch=2.27, r_add=0.160, r_duct=0.825)
