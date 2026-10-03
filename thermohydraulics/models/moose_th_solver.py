"""
MOOSE / CARDINAL THERMAL SOLVER INTERFACE AND PIN BENCHMARK MODULE
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/models/moose_th_solver.py

Provides:
1. High-fidelity reduced thermohydraulic solver for single pins and full core:
   - 1D axial helium enthalpy heating
   - Local axial fluid properties, Re(z), Pr(z), Nu(z), h(z), f(z), dp(z)
   - Nb-1Zr cladding conduction with k(T) = 41.8 + 0.0125 * T
   - Pellet-cladding interface gap conductance model
   - UN fuel conduction with k(T) = 1.37 * T^0.41
   - Exact radial integration for volume-averaged fuel temperature
2. Cardinal / MOOSE finite-element qualification module:
   - Exports 2D RZ axisymmetric MOOSE input files
   - Executes Cardinal (cardinal-opt) for conjugate heat transfer
   - Extracts postprocessors and compares with reduced model.
"""

import os
import sys
import math
import subprocess
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/aleri/OpenMC_projects/first_test/thermohydraulics")
from models.th_properties import (
    R_FUEL, T_CLAD, R_CLAD_OUT, T_GAP, R_CHAN_OUT, H_CORE,
    D_H, A_FLOW_PIN, Q_TOTAL, N_PINS, T_INLET, DELTA_T_BULK,
    cp_helium, pin_coolant_mass_flow, rho_helium, h_convection,
    k_fuel_un, k_clad_nb1zr, reynolds_helium, prandtl_helium,
    nusselt_annulus, friction_factor_annulus, H_GAP_NOMINAL
)

CARDINAL_EXEC = "/home/aleri/projects/cardinal/cardinal-opt"


def solve_pin_thermal_profile(q_vol_slices,
                              m_dot_pin,
                              z_edges,
                              h_gap=H_GAP_NOMINAL,
                              htc_multiplier=1.0,
                              htc_model="gnielinski_trans"):
    """
    Solves axial thermohydraulics and radial heat conduction for a fuel pin across N_z axial slices.
    
    Returns a dictionary of detailed local fields:
    - T_bulk: bulk helium temperature [K]
    - T_clad_out: outer cladding surface temperature [K]
    - T_clad_in: inner cladding surface temperature [K]
    - T_fuel_surf: fuel pellet outer surface temperature [K]
    - T_fuel_center: fuel centerline temperature [K]
    - T_fuel_vol: volume-averaged fuel temperature per slice [K] (for Doppler feedback)
    - h_conv: heat transfer coefficient array [W / (m2 K)]
    - Re: local Reynolds number array
    - Nu: local Nusselt number array
    - Pr: local Prandtl number array
    - rho_cool: local helium density [kg / m3]
    - dp_fric_slice: frictional pressure drop across slice [Pa]
    - dp_accel_slice: acceleration pressure drop across slice [Pa]
    """
    n_slices = len(q_vol_slices)
    dz = H_CORE / n_slices
    v_slice = math.pi * (R_FUEL**2) * dz
    a_clad_slice = 2.0 * math.pi * R_CLAD_OUT * dz
    cp = cp_helium()
    G = m_dot_pin / A_FLOW_PIN  # Mass flux [kg / (m2 s)]
    
    T_bulk = np.zeros(n_slices)
    T_clad_out = np.zeros(n_slices)
    T_clad_in = np.zeros(n_slices)
    T_fuel_surf = np.zeros(n_slices)
    T_fuel_center = np.zeros(n_slices)
    T_fuel_vol = np.zeros(n_slices)
    h_conv = np.zeros(n_slices)
    Re_arr = np.zeros(n_slices)
    Nu_arr = np.zeros(n_slices)
    Pr_arr = np.zeros(n_slices)
    rho_cool = np.zeros(n_slices)
    dp_fric_slice = np.zeros(n_slices)
    dp_accel_slice = np.zeros(n_slices)
    
    T_current = T_INLET
    for s in range(n_slices):
        q_vol = q_vol_slices[s]
        p_slice = q_vol * v_slice
        
        # 1. Enthalpy rise across slice
        delta_T_b = p_slice / (m_dot_pin * cp)
        T_bulk_exit = T_current + delta_T_b
        T_bulk_mid = 0.5 * (T_current + T_bulk_exit)
        T_bulk[s] = T_bulk_exit
        
        # 2. Local fluid properties
        rho_s = rho_helium(T_bulk_mid)
        rho_cool[s] = rho_s
        Re_s = reynolds_helium(m_dot_pin, T_bulk_mid)
        Re_arr[s] = Re_s
        Pr_s = prandtl_helium(T_bulk_mid)
        Pr_arr[s] = Pr_s
        Nu_s = nusselt_annulus(Re_s, Pr_s, model=htc_model)
        Nu_arr[s] = Nu_s
        
        # 3. Heat transfer coefficient
        h_s = h_convection(m_dot_pin, T_bulk_mid, htc_multiplier=htc_multiplier, model=htc_model)
        h_conv[s] = h_s
        
        # 4. Pressure drop in slice
        f_s = friction_factor_annulus(Re_s)
        dp_fric_slice[s] = f_s * (dz / D_H) * (G**2) / (2.0 * rho_s)
        rho_in_s = rho_helium(T_current)
        rho_out_s = rho_helium(T_bulk_exit)
        dp_accel_slice[s] = (G**2) * (1.0 / rho_out_s - 1.0 / rho_in_s)
        
        T_current = T_bulk_exit
        
        # 5. Cladding outer surface temperature
        q_flux_out = p_slice / a_clad_slice  # W / m2
        T_co = T_bulk_mid + q_flux_out / h_s
        T_clad_out[s] = T_co
        
        # 6. Cladding conduction with temperature-dependent k_c(T)
        k_c_guess = k_clad_nb1zr(T_co)
        delta_T_c = (p_slice / (2.0 * math.pi * dz)) * math.log(R_CLAD_OUT / R_FUEL) / k_c_guess
        k_c = k_clad_nb1zr(T_co + 0.5 * delta_T_c)
        delta_T_c = (p_slice / (2.0 * math.pi * dz)) * math.log(R_CLAD_OUT / R_FUEL) / k_c
        T_ci = T_co + delta_T_c
        T_clad_in[s] = T_ci
        
        # 7. Pellet-cladding interface gap
        if h_gap is not None and h_gap > 0.0 and h_gap < 1e8:
            delta_T_gap = (p_slice / (2.0 * math.pi * R_FUEL * dz)) / h_gap
        else:
            delta_T_gap = 0.0
        T_fs = T_ci + delta_T_gap
        T_fuel_surf[s] = T_fs
        
        # 8. Fuel pellet conduction with k_f(T) = 1.37 * T^0.41
        # int_Tfs^T0 k_f(T) dT = q''' * r_f^2 / 4 = (p_slice / dz) / (4 * pi)
        integral_k = (p_slice / dz) / (4.0 * math.pi)
        T0_val = (T_fs**1.41 + (1.41 / 1.37) * integral_k)**(1.0 / 1.41)
        T_fuel_center[s] = T0_val
        
        # 9. Volume-averaged fuel temperature via radial integration
        r_arr = np.linspace(0.0, R_FUEL, 51)
        T_r = (T0_val**1.41 - (1.41 / 1.37) * (q_vol * (r_arr**2) / 4.0))**(1.0 / 1.41)
        T_vol = np.trapezoid(2.0 * r_arr * T_r, r_arr) / (R_FUEL**2)
        T_fuel_vol[s] = T_vol
        
    return {
        "T_bulk": T_bulk,
        "T_clad_out": T_clad_out,
        "T_clad_in": T_clad_in,
        "T_fuel_surf": T_fuel_surf,
        "T_fuel_center": T_fuel_center,
        "T_fuel_vol": T_fuel_vol,
        "h_conv": h_conv,
        "Re": Re_arr,
        "Nu": Nu_arr,
        "Pr": Pr_arr,
        "rho_cool": rho_cool,
        "dp_fric_total_Pa": np.sum(dp_fric_slice),
        "dp_accel_total_Pa": np.sum(dp_accel_slice)
    }


def generate_moose_pin_input(pin_name, z_edges, q_vol_slices, h_coeff, T_bulk_slices, output_file):
    """
    Generates a 2D axisymmetric RZ MOOSE input file for a fuel pin with discrete axial heating
    and convective heat flux boundary condition to helium coolant.
    """
    n_slices = len(q_vol_slices)
    dz = H_CORE / n_slices
    
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
            q_x.append(f"{H_CORE:.5f}")
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

    moose_content = f"""# MOOSE 2D Axisymmetric Pin Cell Model: {pin_name}
[Mesh]
  coord_type = RZ
  [fuel_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 60
    xmin = 0.0
    xmax = {R_FUEL:.5f}
    ymin = 0.0
    ymax = {H_CORE:.5f}
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
    nx = 6
    ny = 60
    xmin = {R_FUEL:.5f}
    xmax = {R_CLAD_OUT:.5f}
    ymin = 0.0
    ymax = {H_CORE:.5f}
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
  # UN Fuel: k = 1.37 * T^0.41 [Hayes et al. 1990]
  [k_fuel]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = 'T'
    expression = '1.37 * (T^0.41)'
    block = '1'
  []
  # Nb-1Zr Clad: k = 41.8 + 0.0125 * T [Touloukian 1970]
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
  [Q_removed_boundary]
    type = SideDiffusiveFluxIntegral
    diffusivity = thermal_conductivity
    variable = T
    boundary = 'clad_right'
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
  [T_fuel_exit_center]
    type = PointValue
    variable = T
    point = '0.0 {H_CORE:.5f} 0.0'
  []
  [T_clad_exit_surf]
    type = PointValue
    variable = T
    point = '{R_CLAD_OUT:.5f} {H_CORE:.5f} 0.0'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
[]

[Outputs]
  csv = true
[]
"""
    with open(output_file, "w") as f:
        f.write(moose_content)


def run_moose_pin(input_file):
    """Executes Cardinal on a MOOSE input file and returns postprocessor dictionary."""
    cmd = [CARDINAL_EXEC, "-i", input_file]
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Cardinal execution failed on {input_file}:\n{result.stderr}")
    
    csv_file = input_file.replace(".i", "_out.csv")
    df = pd.read_csv(csv_file)
    last_row = df.iloc[-1].to_dict()
    return last_row
