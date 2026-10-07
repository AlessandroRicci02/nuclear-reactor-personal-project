#!/usr/bin/env python3
"""BOL screening model for third_try: TB35 lattice with a Be-metal reflector.

The BOL result is a neutronic filter only. The TB35 hot/cold thermal architecture
has not been revised or shown to satisfy a 1150 K outlet and YH2 below 950 K.
The source script was copied from GEOM_OPT15_TB35; some legacy names are retained.
The screened variant uses 217 pins, pitch 2.30 cm, active height 47.0 cm,
radial Be thickness 8.8 cm, axial Be thickness 8.5 cm, and ten B4C drums.
Temperatures used in the fixed BOL model are assumptions, not thermal results.
"""

import os
import sys
import math
import numpy as np
import openmc

OPENMC_EXEC = "/home/aleri/miniforge3/envs/openmc-env/bin/openmc"
CROSS_SECTIONS_PATH = "/home/aleri/OpenMC_projects/first_test/new_possible_best_configurations_failed/campaign/data/endfb-viii.0-hdf5/cross_sections.xml"

# Material densities [g/cm3]
RHO_UN = 13.50
RHO_CLAD = 8.57
RHO_YH2 = 4.30
RHO_BEO = 1.85  # Be metal; legacy identifier retained for mass-code compatibility.
RHO_B4C = 2.52
RHO_AL2O3 = 3.95

# Thermophysical constants
P_NOM = 2.50e6          # 2.5 MPa
R_GAS = 8314.462618
M_HE = 4.002602
N_PINS = 217
DRUM_COUNT = 10

# Geometry Constants (GEOM_OPT15_TB35)
H_CORE = 45.50          # cm
R_FUEL = 0.585          # cm
T_CLAD = 0.048          # cm
T_GAP = 0.084           # cm
T_AL2O3 = 0.010         # cm (100 um)
PITCH = 2.280           # cm
REFL_AX_TH = 8.80       # cm
REFL_RAD_TH = 9.30      # cm
DRUM_R = 3.95           # cm
ABS_TH = 1.50           # cm
ABS_ANGLE = 148.0       # deg
DRUM_CLEARANCE = 0.180  # cm
R_ADD_NOM = 0.150       # cm (1.50 mm)

# Nominal operating temperatures [K]
T_INLET_NOM = 800.0
T_OUTLET_NOM = 1150.0
T_MOD_NOM = 875.0
T_FUEL_NOM = 1300.0
T_CLAD_NOM = 1200.0
T_COOL_NOM = 975.0
T_REFL_NOM = 850.0


def rho_helium(P, T):
    return (P * M_HE) / (R_GAS * max(T, 100.0))


def compute_mass_breakdown(r_add=R_ADD_NOM):
    """
    Computes exact reactor subassembly mass breakdown for GEOM_OPT15_TB35.
    """
    r_f = R_FUEL
    r_c = r_f + T_CLAD
    r_ch = r_c + T_GAP
    r_al = r_ch + T_AL2O3

    # UN fuel
    V_fuel = N_PINS * math.pi * r_f**2 * H_CORE
    M_un = V_fuel * RHO_UN / 1000.0
    M_hm = M_un * (238.0 / (238.0 + 14.0))
    u235_kg = M_hm * (19 * 19.60 + 198 * 19.75) / (217 * 100.0)

    # Cladding Nb-1Zr
    V_clad = N_PINS * math.pi * (r_c**2 - r_f**2) * H_CORE
    M_clad = V_clad * RHO_CLAD / 1000.0

    # Al2O3 Liner
    V_al = N_PINS * math.pi * (r_al**2 - r_ch**2) * H_CORE
    M_al = V_al * RHO_AL2O3 / 1000.0

    # Moderator YH2
    A_hex = (math.sqrt(3)/2.0) * PITCH**2
    A_chan = math.pi * r_al**2
    A_aux = 6 * math.pi * (r_add**2)
    A_mod_cell = max(0.0, A_hex - A_chan - A_aux)
    V_mod = N_PINS * A_mod_cell * H_CORE
    M_yh2_matrix = V_mod * RHO_YH2 / 1000.0

    # Core and reflector radii
    R_core = 19.45 * (PITCH / 2.270)
    R_ref = R_core + REFL_RAD_TH
    V_mod_peripheral = max(0.0, math.pi * R_core**2 - N_PINS * A_hex) * H_CORE
    M_yh2_peripheral = V_mod_peripheral * RHO_YH2 / 1000.0
    M_yh2 = M_yh2_matrix + M_yh2_peripheral

    # Absorber B4C in drums
    A_drum = math.pi * DRUM_R**2
    A_abs = (ABS_ANGLE / 360.0) * math.pi * (DRUM_R**2 - (DRUM_R - ABS_TH)**2)
    A_beo_drum = A_drum - A_abs
    V_abs_tot = DRUM_COUNT * A_abs * H_CORE
    M_b4c = V_abs_tot * RHO_B4C / 1000.0
    V_beo_drums = DRUM_COUNT * A_beo_drum * H_CORE

    # Reflector BeO
    V_rad_annulus = math.pi * (R_ref**2 - R_core**2) * H_CORE - DRUM_COUNT * A_drum * H_CORE + V_beo_drums
    # The axial central duct overlaps each of the six auxiliary channels.
    # Subtract their union rather than subtracting the overlap twice.
    ax_duct_r = 0.860
    aux_offset = 0.910 * (PITCH / 2.270)
    if aux_offset >= ax_duct_r + r_add:
        overlap = 0.0
    elif aux_offset <= abs(ax_duct_r - r_add):
        overlap = math.pi * min(ax_duct_r, r_add)**2
    else:
        d, big_r, small_r = aux_offset, ax_duct_r, r_add
        overlap = (big_r**2 * math.acos((d*d + big_r**2 - small_r**2)/(2*d*big_r))
                   + small_r**2 * math.acos((d*d + small_r**2 - big_r**2)/(2*d*small_r))
                   - .5 * math.sqrt((-d+big_r+small_r)*(d+big_r-small_r)
                                      *(d-big_r+small_r)*(d+big_r+small_r)))
    end_holes_area = N_PINS * (math.pi*ax_duct_r**2 + 6*math.pi*r_add**2 - 6*overlap)
    V_ax = 2.0 * (math.pi*R_ref**2 - end_holes_area) * REFL_AX_TH
    V_beo_tot = V_rad_annulus + V_ax
    M_beo = V_beo_tot * RHO_BEO / 1000.0

    M_tot = M_un + M_clad + M_al + M_yh2 + M_beo + M_b4c

    return {
        "r_add_cm": r_add,
        "r_add_mm": r_add * 10.0,
        "M_un_kg": M_un,
        "M_hm_kg": M_hm,
        "u235_kg": u235_kg,
        "M_clad_kg": M_clad,
        "M_al_kg": M_al,
        "M_yh2_kg": M_yh2,
        "M_yh2_matrix_kg": M_yh2_matrix,
        "M_yh2_peripheral_kg": M_yh2_peripheral,
        "M_beo_kg": M_beo,
        "M_b4c_kg": M_b4c,
        "M_tot_kg": M_tot,
        "R_core_cm": R_core,
        "R_ref_cm": R_ref,
        "margin_to_600kg": 600.0 - M_tot
    }


def build_openmc_model_tb35(run_dir,
                            r_fuel=R_FUEL,
                            clad_th=T_CLAD,
                            chan_gap=T_GAP,
                            al2o3_mod_th=T_AL2O3,
                            pitch=PITCH,
                            core_h=H_CORE,
                            refl_ax_th=REFL_AX_TH,
                            refl_rad_th=REFL_RAD_TH,
                            enr_inner=19.60,
                            enr_outer=19.75,
                            drum_r=DRUM_R,
                            abs_th=ABS_TH,
                            abs_angle=ABS_ANGLE,
                            drum_clearance=DRUM_CLEARANCE,
                            drum_angles=None,
                            reflector_kind="Be",
                            matrix_kind="YH2",
                            peripheral_kind="YH2",
                            r_add=R_ADD_NOM,
                            T_fuel_grid=None,
                            T_clad_grid=None,
                            T_cool_grid=None,
                            T_mod=T_MOD_NOM,
                            T_refl=T_REFL_NOM,
                            T_fuel_offset=0.0,
                            T_mod_offset=0.0,
                            T_liner_offset=0.0,
                            T_clad_offset=0.0,
                            T_refl_offset=0.0,
                            cool_density_mult=1.0,
                            mod_density_mult=1.0,
                            f_in_mat=None,
                            f_out_mat=None,
                            n_axial_slices=10,
                            particles=150000,
                            batches=45,
                            inactive=15,
                            seed=42,
                            tallies_enabled=False):
    """
    Builds the high-fidelity OpenMC model for GEOM_OPT15_TB35.
    """
    os.makedirs(run_dir, exist_ok=True)
    openmc.reset_auto_ids()
    openmc.config['cross_sections'] = CROSS_SECTIONS_PATH

    if drum_angles is None:
        drum_angles = [0.0] * DRUM_COUNT
    if matrix_kind not in ("YH2", "Be") or peripheral_kind not in ("YH2", "Be"):
        raise ValueError("matrix_kind and peripheral_kind must be YH2 or Be")

    # Axial profile for coolant bulk temperature (linear from 892 to 1110 K)
    z_mid_frac = (np.arange(n_axial_slices) + 0.5) / n_axial_slices
    T_cool_axial_nom = T_INLET_NOM + z_mid_frac * (T_OUTLET_NOM - T_INLET_NOM)

    if T_fuel_grid is None:
        # Chopped cosine axial profile centered at 1195 K
        z_axial = np.linspace(-core_h/2.0, core_h/2.0, n_axial_slices)
        cos_shape = np.cos(z_axial * (math.pi / (core_h * 1.35)))
        cos_shape = cos_shape / np.mean(cos_shape)
        T_fuel_grid = np.zeros((N_PINS, n_axial_slices))
        for p in range(N_PINS):
            T_fuel_grid[p, :] = T_FUEL_NOM + 45.0 * (cos_shape - 1.0)

    if T_clad_grid is None:
        T_clad_grid = np.zeros((N_PINS, n_axial_slices))
        for p in range(N_PINS):
            T_clad_grid[p, :] = T_CLAD_NOM + 30.0 * (z_mid_frac - 0.5)

    if T_cool_grid is None:
        T_cool_grid = np.tile(T_cool_axial_nom, (N_PINS, 1))

    T_fuel_grid_eff = T_fuel_grid + T_fuel_offset
    T_clad_grid_eff = T_clad_grid + T_clad_offset
    T_mod_eff = T_mod + T_mod_offset
    T_refl_eff = T_refl + T_refl_offset

    r_f = r_fuel
    r_c = r_f + clad_th
    r_ch = r_c + chan_gap
    r_m_in = r_ch + al2o3_mod_th

    r_add_pos = 0.910 * (pitch / 2.270)

    core_r = 19.45 * (pitch / 2.270)
    refl_rad_outer_r = core_r + refl_rad_th
    total_h = core_h + 2.0 * refl_ax_th
    ax_duct_r = 0.860
    drum_rad_pos = core_r + drum_clearance + drum_r

    # 1. Materials
    v_in = 19.0 * math.pi * (r_f**2) * core_h
    v_out = 198.0 * math.pi * (r_f**2) * core_h

    if f_in_mat is not None:
        f_in = f_in_mat.clone()
        f_in.id = 1
    else:
        f_in = openmc.Material(material_id=1, name="UN Fuel Inner (19.60%)")
        f_in.set_density("g/cm3", RHO_UN)
        f_in.add_nuclide("U235", (enr_inner / 100.0) * (238.0 / (238.0 + 14.0)), "wo")
        f_in.add_nuclide("U238", (1.0 - enr_inner / 100.0) * (238.0 / (238.0 + 14.0)), "wo")
        f_in.add_nuclide("N14", 14.0 / (238.0 + 14.0), "wo")
        f_in.temperature = float(np.mean(T_fuel_grid_eff[:19, :]))
    f_in.volume = v_in
    f_in.depletable = True

    if f_out_mat is not None:
        f_out = f_out_mat.clone()
        f_out.id = 2
    else:
        f_out = openmc.Material(material_id=2, name="UN Fuel Outer (19.75%)")
        f_out.set_density("g/cm3", RHO_UN)
        f_out.add_nuclide("U235", (enr_outer / 100.0) * (238.0 / (238.0 + 14.0)), "wo")
        f_out.add_nuclide("U238", (1.0 - enr_outer / 100.0) * (238.0 / (238.0 + 14.0)), "wo")
        f_out.add_nuclide("N14", 14.0 / (238.0 + 14.0), "wo")
        f_out.temperature = float(np.mean(T_fuel_grid_eff[19:, :]))
    f_out.volume = v_out
    f_out.depletable = True

    clad = openmc.Material(material_id=3, name="Nb-1Zr Cladding")
    clad.set_density("g/cm3", RHO_CLAD)
    clad.add_element("Nb", 0.99, "wo")
    clad.add_element("Zr", 0.01, "wo")
    clad.temperature = float(np.mean(T_clad_grid_eff))

    al2o3 = openmc.Material(material_id=4, name="Al2O3 Liner")
    al2o3.set_density("g/cm3", RHO_AL2O3)
    al2o3.add_element("Al", 2.0, "ao")
    al2o3.add_element("O", 3.0, "ao")
    al2o3.temperature = T_MOD_NOM + T_liner_offset

    # Helium Coolant slices (10 axial slices)
    cool_axial_mats = []
    for s in range(n_axial_slices):
        T_cool_s = float(np.mean(T_cool_grid[:, s]))
        rho_he_si = (rho_helium(P_NOM, T_cool_s) / 1000.0) * cool_density_mult
        m_cool_s = openmc.Material(material_id=20 + s, name=f"He Coolant S{s}")
        m_cool_s.set_density("g/cm3", max(rho_he_si, 1e-7))
        m_cool_s.add_nuclide("He4", 1.0, "ao")
        m_cool_s.temperature = T_cool_s
        cool_axial_mats.append(m_cool_s)

    T_plenum = float(np.mean(T_cool_grid))
    rho_plenum = (rho_helium(P_NOM, T_plenum) / 1000.0) * cool_density_mult
    cool_plenum = openmc.Material(material_id=5, name="He Coolant Plenum")
    cool_plenum.set_density("g/cm3", max(rho_plenum, 1e-7))
    cool_plenum.add_nuclide("He4", 1.0, "ao")
    cool_plenum.temperature = T_plenum

    mod = openmc.Material(material_id=6, name="Canned YH2.0")
    mod.set_density("g/cm3", RHO_YH2 * mod_density_mult)
    mod.add_nuclide("Y89", 1.0, "ao")
    mod.add_nuclide("H1", 2.0, "ao")
    mod.add_s_alpha_beta("c_H_in_YH2")
    mod.add_s_alpha_beta("c_Y_in_YH2")
    mod.temperature = T_mod_eff

    if reflector_kind == "Be":
        refl = openmc.Material(material_id=7, name="Be Reflector")
        refl.set_density("g/cm3", 1.85)
        refl.add_nuclide("Be9", 1.0, "ao")
        refl.add_s_alpha_beta("c_Be")
    elif reflector_kind == "BeO":
        refl = openmc.Material(material_id=7, name="BeO Reflector")
        refl.set_density("g/cm3", 2.90)
        refl.add_nuclide("Be9", 1.0, "ao")
        refl.add_nuclide("O16", 1.0, "ao")
        refl.add_s_alpha_beta("c_Be_in_BeO")
        refl.add_s_alpha_beta("c_O_in_BeO")
    else:
        raise ValueError(f"Unsupported reflector_kind: {reflector_kind}")
    refl.temperature = T_refl_eff

    absorb = openmc.Material(material_id=8, name="B4C Absorber")
    absorb.set_density("g/cm3", RHO_B4C)
    absorb.add_nuclide("B10", 4.0 * 0.90, "ao")
    absorb.add_nuclide("B11", 4.0 * 0.10, "ao")
    absorb.add_element("C", 1.0, "ao")
    absorb.temperature = T_refl_eff

    all_materials = [f_in, f_out, clad, al2o3, cool_plenum, mod, refl, absorb] + cool_axial_mats
    materials = openmc.Materials(all_materials)

    # 2. Geometry
    z_cb = openmc.ZPlane(z0=-core_h / 2.0)
    z_ct = openmc.ZPlane(z0=+core_h / 2.0)
    z_rb = openmc.ZPlane(z0=-total_h / 2.0, boundary_type="vacuum")
    z_rt = openmc.ZPlane(z0=+total_h / 2.0, boundary_type="vacuum")

    z_slice_edges = np.linspace(-core_h / 2.0, +core_h / 2.0, n_axial_slices + 1)
    z_slice_planes = [openmc.ZPlane(z0=z_val) for z_val in z_slice_edges]

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

    fuel_cells_map = {}
    cool_cells_map = {}
    clad_cells_map = {}
    al2o3_cells_map = {}

    def make_pin_unique(pin_idx, f_mat):
        pin_cells = []
        for s_idx in range(n_axial_slices):
            z_low = z_slice_planes[s_idx]
            z_high = z_slice_planes[s_idx + 1]
            T_f_cell = float(T_fuel_grid_eff[pin_idx, s_idx])
            T_c_cell = float(T_clad_grid_eff[pin_idx, s_idx])
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

            c_al = openmc.Cell(name=f"Al2O3_Mod_P{pin_idx}_S{s_idx}", fill=al2o3, region=+cyl_ch & -cyl_m_in & +z_low & -z_high)
            c_al.temperature = T_MOD_NOM + T_liner_offset
            pin_cells.append(c_al)
            al2o3_cells_map[(pin_idx, s_idx)] = c_al

        aux_cells = []
        aux_union = None
        for c in aux_cyls:
            aux_cells.append(openmc.Cell(fill=cool_plenum, region=-c & +z_cb & -z_ct))
            aux_union = -c if aux_union is None else (aux_union | -c)

        mod_reg = +cyl_m_in & +z_cb & -z_ct
        if aux_union is not None:
            mod_reg = mod_reg & ~aux_union
        c_m = openmc.Cell(fill=mod if matrix_kind == "YH2" else refl,
                          region=mod_reg)

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
        c_bb.temperature = T_refl_eff
        c_bt = openmc.Cell(fill=refl, region=+cyl_ax & +z_ct & -z_rt & ~aux_ax_union)
        c_bt.temperature = T_refl_eff

        u = openmc.Universe(cells=pin_cells + [c_m, c_db, c_bb, c_dt, c_bt] + aux_cells + aux_db + aux_dt)
        return u

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
        openmc.Cell(fill=mod if peripheral_kind == "YH2" else refl,
                    region=+z_cb & -z_ct),
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

    tallies = None
    if tallies_enabled:
        fuel_cell_ids = [c.id for c in fuel_cells_map.values()]
        cell_filter = openmc.CellFilter(fuel_cell_ids)
        t_pin_heat = openmc.Tally(name="fuel_cell_heating")
        t_pin_heat.filters = [cell_filter]
        t_pin_heat.scores = ["heating", "fission"]
        tallies = openmc.Tallies([t_pin_heat])

    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = particles
    settings.batches = batches
    settings.inactive = inactive
    settings.temperature = {"multipole": True, "method": "interpolation", "tolerance": 10.0}
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

if __name__ == "__main__":
    m = compute_mass_breakdown()
    print("=" * 70)
    print("GEOM_OPT15_TB35 MASS BREAKDOWN")
    print("=" * 70)
    for k, v in m.items():
        print(f"  {k:<20}: {v}")
