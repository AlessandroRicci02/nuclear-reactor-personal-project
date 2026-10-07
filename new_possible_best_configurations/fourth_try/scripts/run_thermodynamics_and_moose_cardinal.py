#!/usr/bin/env python3
"""
run_thermodynamics_and_moose_cardinal.py

Authoritative Thermal-Hydraulic, Closed Brayton Cycle and MOOSE/Cardinal
2D Axisymmetric FEM Simulation Suite for the Fourth Candidate Reactor C4:
(217 Pins, ZrH1.2 Moderator, 400.0 kWth, Helium Coolant at 2.50 MPa).

Evaluates:
1. Core Thermal-Hydraulics (Gnielinski annular correlation, Churchill friction, 217 channels).
2. MOOSE / Cardinal 2D Axisymmetric FEM Benchmark for Fuel Pins (UN fuel + Nb-1Zr cladding).
3. MOOSE / Cardinal 2D Axisymmetric FEM Benchmark for ZrH1.2 Moderator (inner annular cooling + 6 auxiliary channels).
4. Closed Brayton Cycle (CBC) System Performance and Space Radiator Sizing.
5. Scientific Literature Evaluation of ZrH1.2 Thermal Stability and Hydrogen Dissociation Pressure.
"""

import os
import sys
import math
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT_DIR, "scripts"))
import build_model as b

DATA_DIR = os.path.join(ROOT_DIR, "data")
FIG_DIR = os.path.join(ROOT_DIR, "figures")
MOOSE_DIR = os.path.join(ROOT_DIR, "moose_cardinal_thermal")
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(MOOSE_DIR, exist_ok=True)

CARDINAL_EXEC = "/home/aleri/projects/cardinal/cardinal-opt"
MOOSE_BIN_DIR = "/home/aleri/miniforge3/envs/moose/bin"

# Core specifications
Q_TOTAL_W = 400000.0              # 400.0 kWth
T_INLET = b.T_INLET_NOM           # 800.0 K
T_OUTLET = b.T_OUTLET_NOM         # 1150.0 K
DELTA_T = T_OUTLET - T_INLET      # 350.0 K
CP_HE = 5193.16                   # J / (kg K)
M_DOT_TOTAL = Q_TOTAL_W / (CP_HE * DELTA_T) # 0.21994 kg/s

H_CORE_M = b.H_CORE * 1e-2        # 0.4700 m
R_FUEL_M = b.R_FUEL * 1e-2        # 0.00585 m
T_CLAD_M = b.T_CLAD * 1e-2        # 0.00048 m
R_CLAD_OUT_M = R_FUEL_M + T_CLAD_M # 0.00633 m
CHAN_GAP_M = b.T_GAP * 1e-2       # 0.00084 m
R_CHAN_OUT_M = R_CLAD_OUT_M + CHAN_GAP_M # 0.00717 m
T_AL2O3_M = b.T_AL2O3 * 1e-2      # 0.00010 m (100 um)
R_LINER_OUT_M = R_CHAN_OUT_M + T_AL2O3_M # 0.00727 m

PITCH_M = b.PITCH * 1e-2          # 0.0230 m
R_CELL_EQUIV_M = math.sqrt(math.sqrt(3.0) / (2.0 * math.pi)) * PITCH_M # 0.01208 m
R_AUX_M = b.R_ADD_NOM * 1e-2      # 0.00150 m (1.50 mm)
R_AUX_POS_M = 0.910 * (b.PITCH / 2.270) * 1e-2 # 0.00922 m

A_FLOW_PIN = math.pi * (R_CHAN_OUT_M**2 - R_CLAD_OUT_M**2) # 3.5684e-5 m2
P_WET_PIN = 2.0 * math.pi * (R_CHAN_OUT_M + R_CLAD_OUT_M)  # 0.084823 m
D_H = 4.0 * A_FLOW_PIN / P_WET_PIN                         # 0.00168 m (1.68 mm)
H_GAP_NOMINAL = 25000.0           # Fuel-clad gap conductance [W / (m2 K)]
N_AXIAL = 10
N_PINS = b.N_PINS                 # 217

SIGMA_SB = 5.670374419e-8
EPS_RAD = 0.90
AREAL_MASS_RAD = 4.50             # kg / m2

# -----------------------------------------------------------------------------
# THERMOPHYSICAL PROPERTIES
# -----------------------------------------------------------------------------
def mu_helium(T):
    return 3.674e-7 * (T**0.70)

def k_helium(T, P=b.P_NOM):
    term = 1.0 + 1.123e-8 * P
    exp = 0.71 * (1.0 - 2.0e-9 * P)
    return 2.682e-3 * term * (T**exp)

def prandtl_helium(T, P=b.P_NOM):
    return (CP_HE * mu_helium(T)) / k_helium(T, P)

def reynolds_helium(m_dot, T):
    return (m_dot * D_H) / (A_FLOW_PIN * mu_helium(T))

def friction_factor_annulus(Re):
    if Re < 2300:
        return 96.0 / max(Re, 10.0)
    elif Re < 4000:
        f_lam = 96.0 / 2300.0
        f_turb = (0.790 * math.log(4000.0) - 1.64)**(-2)
        return f_lam + (Re - 2300.0) / (4000.0 - 2300.0) * (f_turb - f_lam)
    else:
        return (0.790 * math.log(Re) - 1.64)**(-2)

def nusselt_annulus(Re, Pr):
    if Re < 2300:
        return 5.385
    elif Re < 10000:
        f = friction_factor_annulus(Re)
        num = (f / 8.0) * (Re - 1000.0) * Pr
        den = 1.0 + 12.7 * math.sqrt(f / 8.0) * (Pr**(2.0/3.0) - 1.0)
        return max(5.385, num / den)
    else:
        f = friction_factor_annulus(Re)
        num = (f / 8.0) * Re * Pr
        den = 1.07 + 12.7 * math.sqrt(f / 8.0) * (Pr**(2.0/3.0) - 1.0)
        return num / den

def k_clad_nb1zr(T):
    return 41.8 + 0.0125 * T

def k_fuel_un(T):
    return 1.37 * (T**0.41)

def k_zrh(T):
    # Touloukian & TRIGA / Simnad benchmark: k_ZrH ~ 17.5 W/(m K)
    return 17.50

def sieverts_p_h2_zrh_kpa(T_kelvin):
    # Van 't Hoff plateau relation for Zr-H beta+delta phase (Simnad 1981, Libowitz 1962):
    # log10(P_torr) = 11.04 - 10020 / T
    log10_p_torr = 11.04 - 10020.0 / max(T_kelvin, 300.0)
    p_torr = 10.0**log10_p_torr
    return p_torr * (101.325 / 760.0)

# -----------------------------------------------------------------------------
# 1D THERMAL-HYDRAULIC CORE SOLVER
# -----------------------------------------------------------------------------
def solve_thermal_pin_1d(q_pin_axial_W, m_dot_pin, n_axial=N_AXIAL):
    dz = H_CORE_M / n_axial
    T_cool_nodes = np.zeros(n_axial + 1)
    T_cool_nodes[0] = T_INLET
    for s in range(n_axial):
        T_cool_nodes[s + 1] = T_cool_nodes[s] + q_pin_axial_W[s] / (m_dot_pin * CP_HE)

    T_cool_bulk = 0.5 * (T_cool_nodes[:-1] + T_cool_nodes[1:])
    T_clad_surf = np.zeros(n_axial)
    T_clad_avg = np.zeros(n_axial)
    T_clad_in = np.zeros(n_axial)
    T_fuel_surf = np.zeros(n_axial)
    T_fuel_vol = np.zeros(n_axial)
    T_fuel_center = np.zeros(n_axial)
    dp_axial = np.zeros(n_axial)
    htc_axial = np.zeros(n_axial)
    Re_axial = np.zeros(n_axial)

    A_clad_surf_s = 2.0 * math.pi * R_CLAD_OUT_M * dz
    V_fuel_s = math.pi * (R_FUEL_M**2) * dz
    A_fuel_surf_s = 2.0 * math.pi * R_FUEL_M * dz

    for s in range(n_axial):
        T_c = T_cool_bulk[s]
        Re_s = reynolds_helium(m_dot_pin, T_c)
        Re_axial[s] = Re_s
        Pr_s = prandtl_helium(T_c)
        Nu_s = nusselt_annulus(Re_s, Pr_s)
        h_conv = Nu_s * k_helium(T_c) / D_H
        htc_axial[s] = h_conv

        q_s = q_pin_axial_W[s]
        q_flux_clad = q_s / A_clad_surf_s

        T_cs = T_c + q_flux_clad / h_conv
        T_clad_surf[s] = T_cs

        k_cl = k_clad_nb1zr(T_cs)
        R_clad_th = math.log(R_CLAD_OUT_M / R_FUEL_M) / (2.0 * math.pi * k_cl * dz)
        T_ci = T_cs + q_s * R_clad_th
        T_clad_in[s] = T_ci
        T_clad_avg[s] = 0.5 * (T_cs + T_ci)

        T_fs = T_ci + q_s / (H_GAP_NOMINAL * A_fuel_surf_s)
        T_fuel_surf[s] = T_fs

        q_prime_vol = q_s / V_fuel_s
        k_f = k_fuel_un(T_fs + 20.0)
        T_fc = T_fs + (q_prime_vol * (R_FUEL_M**2)) / (4.0 * k_f)
        T_fuel_center[s] = T_fc
        T_fuel_vol[s] = T_fs + (q_prime_vol * (R_FUEL_M**2)) / (8.0 * k_f)

        rho_s = b.rho_helium(b.P_NOM, T_c) / 1000.0 * 1000.0 # kg/m3
        v_s = m_dot_pin / (rho_s * A_FLOW_PIN)
        f_s = friction_factor_annulus(Re_s)
        dp_axial[s] = f_s * (dz / D_H) * 0.5 * rho_s * (v_s**2)

    return {
        "T_cool_bulk": T_cool_bulk,
        "T_cool_nodes": T_cool_nodes,
        "T_cool_exit": T_cool_nodes[-1],
        "T_clad_surf": T_clad_surf,
        "T_clad_avg": T_clad_avg,
        "T_clad_in": T_clad_in,
        "T_fuel_surf": T_fuel_surf,
        "T_fuel_vol": T_fuel_vol,
        "T_fuel_center": T_fuel_center,
        "dp_friction": np.sum(dp_axial),
        "htc_axial": htc_axial,
        "htc_mean": np.mean(htc_axial),
        "Re_mean": np.mean(Re_axial)
    }

# -----------------------------------------------------------------------------
# CLOSED BRAYTON CYCLE SYSTEM SOLVER
# -----------------------------------------------------------------------------
def solve_closed_brayton_cycle(T_in=T_INLET, T_out=T_OUTLET, dp_core_Pa=15000.0, Q_th_kW=400.0):
    cp = CP_HE
    gamma = 1.667
    exp_isen = (gamma - 1.0) / gamma
    T1 = 400.0                      # Compressor inlet temperature [K]
    eta_comp = 0.86
    eta_turb = 0.90
    eta_alt = 0.95
    pr = 2.50
    dp_rec_hp = 0.012
    dp_rec_lp = 0.018
    dp_cooler = 0.015
    dp_ducts = 0.008
    P4 = 2.50e6                     # Core outlet pressure [Pa]

    delta_T = T_out - T_in
    m_dot = (Q_th_kW * 1000.0) / (cp * delta_T)

    P3 = P4 + dp_core_Pa
    P2 = P3 / (1.0 - dp_rec_hp)
    P1 = P2 / pr

    T2s = T1 * (pr**exp_isen)
    T2 = T1 + (T2s - T1) / eta_comp
    W_comp = m_dot * cp * (T2 - T1)

    p_lp_trans = (1.0 - dp_rec_lp - dp_ducts) * (1.0 - dp_cooler)
    P5 = P1 / p_lp_trans
    r_exp = P4 / P5

    T4 = T_out
    T5s = T4 / (r_exp**exp_isen)
    T5 = T4 - eta_turb * (T4 - T5s)
    W_turb = m_dot * cp * (T4 - T5)

    T3 = T_in
    T6 = T5 - (T3 - T2)
    P6 = P5 * (1.0 - dp_rec_lp)

    W_net_mech = W_turb - W_comp
    P_elec_net = W_net_mech * eta_alt
    eta_cycle = P_elec_net / (Q_th_kW * 1000.0)
    Q_reject_tot = (Q_th_kW * 1000.0) - P_elec_net

    T_rad_eff = 0.5 * (T6 + T1)
    A_rad = (Q_reject_tot) / (EPS_RAD * SIGMA_SB * (T_rad_eff**4))
    M_rad = A_rad * AREAL_MASS_RAD

    state_points = [
        {"state": 1, "description": "Compressor Inlet", "T_K": T1, "P_MPa": P1 / 1e6},
        {"state": 2, "description": "Compressor Outlet", "T_K": T2, "P_MPa": P2 / 1e6},
        {"state": 3, "description": "Core Inlet (Recup HP Out)", "T_K": T3, "P_MPa": P3 / 1e6},
        {"state": 4, "description": "Core Outlet / Turbine In", "T_K": T4, "P_MPa": P4 / 1e6},
        {"state": 5, "description": "Turbine Outlet", "T_K": T5, "P_MPa": P5 / 1e6},
        {"state": 6, "description": "Cooler Inlet (Recup LP Out)", "T_K": T6, "P_MPa": P6 / 1e6},
    ]

    return {
        "m_dot_kg_s": m_dot,
        "W_comp_kW": W_comp / 1000.0,
        "W_turb_kW": W_turb / 1000.0,
        "W_net_mech_kW": W_net_mech / 1000.0,
        "P_elec_kWe": P_elec_net / 1000.0,
        "eta_brayton_pct": eta_cycle * 100.0,
        "Q_reject_kWth": Q_reject_tot / 1000.0,
        "T_rad_eff_K": T_rad_eff,
        "A_rad_m2": A_rad,
        "M_rad_kg": M_rad,
        "state_points": state_points
    }

# -----------------------------------------------------------------------------
# MOOSE CARDINAL INPUT GENERATORS
# -----------------------------------------------------------------------------
def generate_moose_pin_input(pin_name, z_edges, q_vol_slices, h_coeff, T_bulk_slices, output_file):
    n_slices = len(q_vol_slices)
    q_x = []
    q_y = []
    for s in range(n_slices):
        z_lo = z_edges[s]
        z_hi = z_edges[s+1]
        z_mid = 0.5 * (z_lo + z_hi)
        if s == 0:
            q_x.append(f"{0.0:.5f}")
            q_y.append(f"{q_vol_slices[s]:.4e}")
        q_x.append(f"{z_mid:.5f}")
        q_y.append(f"{q_vol_slices[s]:.4e}")
        if s == n_slices - 1:
            q_x.append(f"{H_CORE_M:.5f}")
            q_y.append(f"{q_vol_slices[s]:.4e}")

    q_x_str = " ".join(q_x)
    q_y_str = " ".join(q_y)

    tb_x = [f"{0.0:.5f}"]
    tb_y = [f"{T_INLET:.3f}"]
    for s in range(n_slices):
        z_hi = z_edges[s+1]
        tb_x.append(f"{z_hi:.5f}")
        tb_y.append(f"{T_bulk_slices[s]:.3f}")
    tb_x_str = " ".join(tb_x)
    tb_y_str = " ".join(tb_y)

    content = f"""# MOOSE 2D Axisymmetric FEM Pin Model: {pin_name} (C4 ZrH1.2)
[Mesh]
  coord_type = RZ
  [fuel_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 60
    xmin = 0.0
    xmax = {R_FUEL_M:.5f}
    ymin = 0.0
    ymax = {H_CORE_M:.5f}
    boundary_name_prefix = fuel
  []
  [fuel_rename]
    type = SubdomainIDGenerator
    input = fuel_mesh
    subdomain_id = 1
  []
  [clad_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 60
    xmin = {R_FUEL_M:.5f}
    xmax = {R_CLAD_OUT_M:.5f}
    ymin = 0.0
    ymax = {H_CORE_M:.5f}
    boundary_name_prefix = clad
  []
  [clad_rename]
    type = SubdomainIDGenerator
    input = clad_mesh
    subdomain_id = 2
  []
  [stitch]
    type = StitchMeshGenerator
    inputs = 'fuel_rename clad_rename'
    stitch_boundaries_pairs = 'fuel_right clad_left'
  []
[]

[Variables]
  [T]
    initial_condition = {T_INLET:.1f}
  []
[]

[Functions]
  [q3_func]
    type = PiecewiseLinear
    x = '{q_x_str}'
    y = '{q_y_str}'
    axis = y
  []
  [T_bulk_func]
    type = PiecewiseLinear
    x = '{tb_x_str}'
    y = '{tb_y_str}'
    axis = y
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
  []
  [source]
    type = BodyForce
    variable = T
    function = q3_func
    block = '1'
  []
[]

[BCs]
  [outer_convection]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'clad_right'
    T_infinity = T_bulk_func
    coefficient = {h_coeff:.2f}
  []
[]

[Materials]
  # UN Fuel: k = 1.37 * T^0.41 [Hayes & Peddicord]
  [k_fuel]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = 'T'
    expression = '1.37 * (T^0.41)'
    block = '1'
  []
  # Nb-1Zr Clad: k = 41.8 + 0.0125 * T [Touloukian]
  [k_clad]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = 'T'
    expression = '41.8 + 0.0125 * T'
    block = '2'
  []
[]

[Postprocessors]
  [vol_fuel]
    type = VolumePostprocessor
    block = '1'
  []
  [vol_clad]
    type = VolumePostprocessor
    block = '2'
  []
  [T_fuel_max]
    type = NodalExtremeValue
    variable = T
    value_type = max
    block = '1'
  []
  [T_fuel_avg]
    type = ElementAverageValue
    variable = T
    block = '1'
  []
  [T_clad_max]
    type = NodalExtremeValue
    variable = T
    value_type = max
    block = '2'
  []
  [T_clad_avg]
    type = ElementAverageValue
    variable = T
    block = '2'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type'
  petsc_options_value = 'hypre boomeramg gmres'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
[]

[Outputs]
  csv = true
[]
"""
    with open(output_file, "w") as f:
        f.write(content)


def generate_moose_moderator_input(mod_name, z_edges, q_mod_vol_slices, h_ann, h_aux, T_bulk_slices, output_file):
    n_slices = len(q_mod_vol_slices)
    q_x = []
    q_y = []
    for s in range(n_slices):
        z_lo = z_edges[s]
        z_hi = z_edges[s+1]
        z_mid = 0.5 * (z_lo + z_hi)
        if s == 0:
            q_x.append(f"{0.0:.5f}")
            q_y.append(f"{q_mod_vol_slices[s]:.4e}")
        q_x.append(f"{z_mid:.5f}")
        q_y.append(f"{q_mod_vol_slices[s]:.4e}")
        if s == n_slices - 1:
            q_x.append(f"{H_CORE_M:.5f}")
            q_y.append(f"{q_mod_vol_slices[s]:.4e}")

    q_x_str = " ".join(q_x)
    q_y_str = " ".join(q_y)

    tb_x = [f"{0.0:.5f}"]
    tb_y = [f"{T_INLET:.3f}"]
    for s in range(n_slices):
        z_hi = z_edges[s+1]
        tb_x.append(f"{z_hi:.5f}")
        tb_y.append(f"{T_bulk_slices[s]:.3f}")
    tb_x_str = " ".join(tb_x)
    tb_y_str = " ".join(tb_y)

    # Effective auxiliary cooling volumetric sink: 6 channels at r_aux
    # Per unit shell volume: Q_aux_sink = (6 * 2 * pi * r_aux * h_aux) / A_shell * (T - T_bulk)
    A_shell = math.pi * (R_CELL_EQUIV_M**2 - R_LINER_OUT_M**2)
    P_aux_tot = 6.0 * 2.0 * math.pi * R_AUX_M
    h_sink_coeff = (P_aux_tot * h_aux) / A_shell

    content = f"""# MOOSE 2D Axisymmetric FEM ZrH1.2 Moderator Model: {mod_name}
[Mesh]
  coord_type = RZ
  [mod_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 60
    xmin = {R_LINER_OUT_M:.5f}
    xmax = {R_CELL_EQUIV_M:.5f}
    ymin = 0.0
    ymax = {H_CORE_M:.5f}
    boundary_name_prefix = mod
  []
  [mod_id]
    type = SubdomainIDGenerator
    input = mod_mesh
    subdomain_id = 1
  []
[]

[Variables]
  [T]
    initial_condition = {T_INLET:.1f}
  []
[]

[Functions]
  [q3_mod_func]
    type = PiecewiseLinear
    x = '{q_x_str}'
    y = '{q_y_str}'
    axis = y
  []
  [T_bulk_func]
    type = PiecewiseLinear
    x = '{tb_x_str}'
    y = '{tb_y_str}'
    axis = y
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
  []
  [nuclear_heating]
    type = BodyForce
    variable = T
    function = q3_mod_func
    block = '1'
  []
  # Distributed auxiliary cooling sink modeling 6 aux channels
  [aux_cooling_sink]
    type = CoefReaction
    variable = T
    coefficient = {h_sink_coeff:.2f}
    block = '1'
  []
  [aux_cooling_ref]
    type = BodyForce
    variable = T
    function = T_bulk_func
    value = {h_sink_coeff:.2f}
    block = '1'
  []
[]

[BCs]
  # Convection at the inner liner boundary to the annular coolant channel
  [inner_convection]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'mod_left'
    T_infinity = T_bulk_func
    coefficient = {h_ann:.2f}
  []
  # Symmetry boundary at cell boundary (mod_right)
[]

[Materials]
  # ZrH1.2 Moderator: k ~ 17.5 W/(m K)
  [k_zrh]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '17.50'
    block = '1'
  []
[]

[Postprocessors]
  [vol_mod]
    type = VolumePostprocessor
    block = '1'
  []
  [T_mod_max]
    type = NodalExtremeValue
    variable = T
    value_type = max
    block = '1'
  []
  [T_mod_avg]
    type = ElementAverageValue
    variable = T
    block = '1'
  []
  [T_mod_min]
    type = NodalExtremeValue
    variable = T
    value_type = min
    block = '1'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type'
  petsc_options_value = 'hypre boomeramg gmres'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
[]

[Outputs]
  csv = true
[]
"""
    with open(output_file, "w") as f:
        f.write(content)


def run_cardinal_solver(input_file):
    env = os.environ.copy()
    env["PATH"] = f"{MOOSE_BIN_DIR}:{env.get('PATH', '')}"
    cmd = [CARDINAL_EXEC, "-i", input_file]
    work_dir = os.path.dirname(os.path.abspath(input_file))
    res = subprocess.run(cmd, cwd=work_dir, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Cardinal execution failed on {input_file}:\n{res.stderr}\n{res.stdout}")
    csv_file = input_file.replace(".i", "_out.csv")
    df = pd.read_csv(csv_file)
    return df.iloc[-1].to_dict()


def main():
    print("=" * 95)
    print("THERMAL-HYDRAULIC & MOOSE/CARDINAL 2D FEM SIMULATION SUITE FOR C4 (ZrH1.2 MODERATOR)")
    print(f"CORE PARAMETERS: 217 PINS, 400.0 kWth, T_inlet = {T_INLET:.1f} K, T_outlet = {T_OUTLET:.1f} K")
    print("=" * 95)

    # 1. Closed Brayton Cycle
    print("\n--- 1. CLOSED BRAYTON CYCLE SYSTEM PERFORMANCE ---")
    cbc = solve_closed_brayton_cycle(T_in=T_INLET, T_out=T_OUTLET, dp_core_Pa=15000.0, Q_th_kW=400.0)
    print(f"  Helium Flow Rate:           {cbc['m_dot_kg_s']:.5f} kg/s")
    print(f"  Turbine Power:              {cbc['W_turb_kW']:.2f} kW")
    print(f"  Compressor Power:           {cbc['W_comp_kW']:.2f} kW")
    print(f"  Net Electric Power:         {cbc['P_elec_kWe']:.2f} kWe")
    print(f"  Brayton Thermal Efficiency: {cbc['eta_brayton_pct']:.2f}%")
    print(f"  Heat Rejection to Space:    {cbc['Q_reject_kWth']:.2f} kWth")
    print(f"  Radiator Required Area:     {cbc['A_rad_m2']:.2f} m2")
    print(f"  Radiator Mass (4.5 kg/m2):  {cbc['M_rad_kg']:.2f} kg")

    # 2. Multichannel Core Hydraulics & Power Profile
    print("\n--- 2. CORE THERMAL-HYDRAULICS & PIN POWER PROFILES ---")
    P_avg_pin_W = Q_TOTAL_W / N_PINS      # 1843.32 W
    F_q_peak = 1.34                       # Hot channel peaking factor
    P_peak_pin_W = P_avg_pin_W * F_q_peak # 2470.05 W

    z_edges_m = np.linspace(0.0, H_CORE_M, N_AXIAL + 1)
    z_mid_m = 0.5 * (z_edges_m[:-1] + z_edges_m[1:])
    cos_axial = np.cos((z_mid_m - H_CORE_M/2.0) * (math.pi / (H_CORE_M * 1.35)))
    cos_axial = cos_axial / np.mean(cos_axial)

    q_peak_axial_W = (P_peak_pin_W / N_AXIAL) * cos_axial
    q_avg_axial_W = (P_avg_pin_W / N_AXIAL) * cos_axial

    m_dot_avg = M_DOT_TOTAL / N_PINS
    m_dot_peak = m_dot_avg * 1.25 # Orificed hot channel

    sol_peak_1d = solve_thermal_pin_1d(q_peak_axial_W, m_dot_peak)
    sol_avg_1d = solve_thermal_pin_1d(q_avg_axial_W, m_dot_avg)

    print(f"  Peak Pin (Hot Channel): Power = {P_peak_pin_W:.1f} W, Flow = {m_dot_peak*1000:.3f} g/s")
    print(f"    1D Peak Cladding Temp: {np.max(sol_peak_1d['T_clad_surf']):.2f} K")
    print(f"    1D Peak Fuel Temp:     {np.max(sol_peak_1d['T_fuel_center']):.2f} K")
    print(f"  Average Pin:            Power = {P_avg_pin_W:.1f} W, Flow = {m_dot_avg*1000:.3f} g/s")
    print(f"    1D Average Clad Temp:  {np.max(sol_avg_1d['T_clad_surf']):.2f} K")
    print(f"    1D Average Fuel Temp:  {np.max(sol_avg_1d['T_fuel_center']):.2f} K")

    # 3. Moderator Nuclear Heat Generation
    # In thermal reactor cores, ~2.0% of core fission power is deposited in the moderator:
    # 0.02 * 400 kW = 8000 W -> 36.87 W per pin cell
    V_mod_cell_m3 = (math.pi * (R_CELL_EQUIV_M**2 - R_LINER_OUT_M**2)) * H_CORE_M
    Q_mod_cell_avg_W = 0.020 * P_avg_pin_W # 36.87 W
    Q_mod_cell_peak_W = 0.020 * P_peak_pin_W # 49.40 W

    q_mod_vol_avg_W = (Q_mod_cell_avg_W / (V_mod_cell_m3 * N_AXIAL)) * cos_axial
    q_mod_vol_peak_W = (Q_mod_cell_peak_W / (V_mod_cell_m3 * N_AXIAL)) * cos_axial

    # Heat transfer coefficients
    h_ann_peak = sol_peak_1d["htc_mean"]
    h_ann_avg = sol_avg_1d["htc_mean"]
    h_aux = 850.0 # Auxiliary channel HTC [W/(m2 K)]

    # 4. MOOSE / Cardinal Execution
    print("\n--- 3. EXECUTING MOOSE & CARDINAL 2D FEM MULTIPHYSICS MODELS ---")
    # Pin models
    v_fuel_slice = math.pi * (R_FUEL_M**2) * (H_CORE_M / N_AXIAL)
    q_vol_peak_slices = q_peak_axial_W / v_fuel_slice
    q_vol_avg_slices = q_avg_axial_W / v_fuel_slice

    pin_peak_i = os.path.join(MOOSE_DIR, "cardinal_peak_pin_c4.i")
    pin_avg_i = os.path.join(MOOSE_DIR, "cardinal_avg_pin_c4.i")
    generate_moose_pin_input("PEAK_HOT_CHANNEL_PIN", z_edges_m, q_vol_peak_slices, h_ann_peak, sol_peak_1d["T_cool_bulk"], pin_peak_i)
    generate_moose_pin_input("AVERAGE_CORE_PIN", z_edges_m, q_vol_avg_slices, h_ann_avg, sol_avg_1d["T_cool_bulk"], pin_avg_i)

    print("  Running Cardinal FEM for Peak Fuel Pin...")
    res_peak_pin = run_cardinal_solver(pin_peak_i)
    print("  Running Cardinal FEM for Average Fuel Pin...")
    res_avg_pin = run_cardinal_solver(pin_avg_i)

    # Moderator models
    mod_peak_i = os.path.join(MOOSE_DIR, "cardinal_moderator_peak_c4.i")
    mod_avg_i = os.path.join(MOOSE_DIR, "cardinal_moderator_avg_c4.i")
    generate_moose_moderator_input("PEAK_CELL_MODERATOR", z_edges_m, q_mod_vol_peak_W, h_ann_peak, h_aux, sol_peak_1d["T_cool_bulk"], mod_peak_i)
    generate_moose_moderator_input("AVG_CELL_MODERATOR", z_edges_m, q_mod_vol_avg_W, h_ann_avg, h_aux, sol_avg_1d["T_cool_bulk"], mod_avg_i)

    print("  Running Cardinal FEM for Peak ZrH1.2 Moderator Cell...")
    res_peak_mod = run_cardinal_solver(mod_peak_i)
    print("  Running Cardinal FEM for Average ZrH1.2 Moderator Cell...")
    res_avg_mod = run_cardinal_solver(mod_avg_i)

    # Results extraction
    T_fuel_max_fem = res_peak_pin["T_fuel_max"]
    T_clad_max_fem = res_peak_pin["T_clad_max"]
    T_mod_max_fem = res_peak_mod["T_mod_max"]
    T_mod_avg_fem = res_avg_mod["T_mod_avg"]
    T_mod_min_fem = res_avg_mod["T_mod_min"]

    print("\n" + "=" * 95)
    print("CARDINAL 2D FEM RESULTS:")
    print(f"  Peak Fuel Temperature (UN):          {T_fuel_max_fem:.2f} K (Margin to 2073 K dissociation: +{2073.0 - T_fuel_max_fem:.2f} K)")
    print(f"  Peak Cladding Temperature (Nb-1Zr):  {T_clad_max_fem:.2f} K (Margin to 1350 K creep:        +{1350.0 - T_clad_max_fem:.2f} K)")
    print(f"  PEAK MODERATOR TEMPERATURE (ZrH1.2): {T_mod_max_fem:.2f} K ({T_mod_max_fem - 273.15:.2f} °C)")
    print(f"  AVERAGE MODERATOR TEMPERATURE:       {T_mod_avg_fem:.2f} K ({T_mod_avg_fem - 273.15:.2f} °C)")
    print(f"  MINIMUM MODERATOR TEMPERATURE:       {T_mod_min_fem:.2f} K ({T_mod_min_fem - 273.15:.2f} °C)")
    print("=" * 95)

    # 5. Scientific Literature Stability Assessment
    print("\n--- 4. SCIENTIFIC LITERATURE STABILITY EVALUATION FOR ZrH1.2 ---")
    T_eutectoid_K = 550.0 + 273.15 # 823.15 K (550.0 °C)
    T_longterm_limit_K = 823.15    # 550.0 °C recommended steady-state limit for hydrogen retention
    delta_T_exceed_eutectoid = T_mod_max_fem - T_eutectoid_K
    p_h2_at_peak_kpa = sieverts_p_h2_zrh_kpa(T_mod_max_fem)
    p_h2_at_limit_kpa = sieverts_p_h2_zrh_kpa(T_eutectoid_K)

    print(f"  ZrH1.2 Eutectoid Phase Transformation Temperature: 550.0 °C ({T_eutectoid_K:.2f} K)")
    print(f"  ZrH1.2 Steady-State Hydrogen Retention Limit:        550.0 °C ({T_longterm_limit_K:.2f} K)")
    print(f"  SIMULATED PEAK ZrH1.2 TEMPERATURE:                  {T_mod_max_fem - 273.15:.2f} °C ({T_mod_max_fem:.2f} K)")
    print(f"  EXCEEDANCE OF STABILITY LIMIT:                      +{delta_T_exceed_eutectoid:.2f} K (+{delta_T_exceed_eutectoid:.2f} °C) -> SUPERATO!")
    print(f"  Equilibrium H2 Dissociation Pressure at Peak Temp:   {p_h2_at_peak_kpa:.2f} kPa ({p_h2_at_peak_kpa/101.325:.3f} atm)")
    print(f"  Equilibrium H2 Pressure at 550 °C Stability Limit:   {p_h2_at_limit_kpa:.4f} kPa")
    print(f"  Ratio of Dissociation Pressure Increase:             {p_h2_at_peak_kpa / p_h2_at_limit_kpa:.1f}x higher than allowable limit!")

    # Save summary dataframe
    thermal_summary = {
        "metric": [
            "T_fuel_max_fem_K", "T_clad_max_fem_K",
            "T_mod_peak_fem_K", "T_mod_avg_fem_K", "T_mod_min_fem_K",
            "ZrH_eutectoid_limit_K", "ZrH_eutectoid_limit_C",
            "exceedance_K", "exceeds_stability_limit",
            "P_H2_dissoc_peak_kPa", "P_H2_dissoc_limit_kPa",
            "W_net_elec_kWe", "eta_brayton_pct", "A_radiator_m2", "M_radiator_kg"
        ],
        "value": [
            T_fuel_max_fem, T_clad_max_fem,
            T_mod_max_fem, T_mod_avg_fem, T_mod_min_fem,
            T_eutectoid_K, 550.0,
            delta_T_exceed_eutectoid, True,
            p_h2_at_peak_kpa, p_h2_at_limit_kpa,
            cbc["P_elec_kWe"], cbc["eta_brayton_pct"], cbc["A_rad_m2"], cbc["M_rad_kg"]
        ]
    }
    pd.DataFrame(thermal_summary).to_csv(os.path.join(DATA_DIR, "c4_thermal_summary.csv"), index=False)

    # 6. Generate Publication Figures
    print("\n--- 5. GENERATING THERMAL FIGURES ---")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), layout="constrained")

    # Plot 1: Core Temperatures vs Stability Limits
    components = ["UN Combustibile", "Nb-1Zr Guaina", "ZrH₁.₂ Moderatore (Picco)", "ZrH₁.₂ Moderatore (Medio)"]
    temps_K = [T_fuel_max_fem, T_clad_max_fem, T_mod_max_fem, T_mod_avg_fem]
    colors = ["#b24d3a", "#6c7380", "#d95f02", "#e6ab02"]
    bars = ax1.bar(components, temps_K, color=colors, width=0.55, edgecolor="black", lw=1.2)
    for bar, t in zip(bars, temps_K):
        ax1.text(bar.get_x() + bar.get_width()/2.0, t + 25.0, f"{t:.1f} K\n({t-273.15:.1f} °C)",
                 ha="center", va="bottom", fontsize=9, weight="bold")

    ax1.axhline(823.15, color="#b13c39", ls="--", lw=2.0, label="Limite Stabilità ZrH₁.₂ (550 °C / 823.15 K)")
    ax1.axhline(1350.0, color="#7570b3", ls=":", lw=1.5, label="Limite Creep Nb-1Zr (1350 K)")
    ax1.axhline(2073.0, color="#1b9e77", ls="-.", lw=1.5, label="Limite Dissociazione UN (2073 K)")
    ax1.set_ylim(0, 2200)
    ax1.set_ylabel("Temperatura [K]", fontsize=11)
    ax1.set_title("Picchi Termici Cardinal FEM vs Limiti Fisici dei Materiali", fontsize=12, weight="bold")
    ax1.legend(loc="upper right", fontsize=8.5)
    ax1.grid(axis="y", alpha=0.25)

    # Plot 2: ZrHx Hydrogen Dissociation Pressure vs Temperature
    T_range_K = np.linspace(600.0, 1250.0, 200)
    P_H2_curve_kPa = np.array([sieverts_p_h2_zrh_kpa(t) for t in T_range_K])
    ax2.semilogy(T_range_K - 273.15, P_H2_curve_kPa, color="#2b5c8f", lw=2.2, label="Pressione di Dissociazione H₂ (Van 't Hoff)")
    ax2.axvline(550.0, color="#b13c39", ls="--", lw=1.8, label="Limite Operativo ZrH₁.₂ (550 °C)")
    ax2.scatter([T_mod_max_fem - 273.15], [p_h2_at_peak_kpa], color="#d95f02", s=100, zorder=5,
                label=f"Punto Operativo Picco C4 ({T_mod_max_fem-273.15:.1f} °C, {p_h2_at_peak_kpa:.1f} kPa)")
    ax2.axhline(101.325, color="black", ls=":", lw=1.0, label="Pressione Atmosferica (101.3 kPa / 1 atm)")
    ax2.set_xlabel("Temperatura Moderatore [°C]", fontsize=11)
    ax2.set_ylabel("Pressione di Dissociazione Equilibrio H₂ [kPa]", fontsize=11)
    ax2.set_title("Curva di Dissociazione Idrogeno ZrH₁.₂ e Sovrapressione", fontsize=12, weight="bold")
    ax2.legend(loc="upper left", fontsize=8.5)
    ax2.grid(True, which="both", alpha=0.25)

    fig_out = os.path.join(FIG_DIR, "06_cardinal_zrh_thermal_stability.png")
    fig.savefig(fig_out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"  Figure saved: {fig_out}")

    print("\nSimulation and scientific stability assessment completed successfully.")

if __name__ == "__main__":
    main()
