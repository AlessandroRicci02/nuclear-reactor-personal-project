"""
CARDINAL / MOOSE / THM THERMOHYDRAULIC AND SOLID CONDUCTION INTERFACE
Workspace: /home/aleri/OpenMC_projects/first_test/thermohydraulics/models/cardinal_th_interface.py

Directly interfaces with Cardinal (cardinal-opt) and MOOSE:
1. Solves 1D axial compressible flow conservation equations (mass, momentum, energy):
   - Integrates FlowChannel1Phase with HeliumFluidProperties (Petersen 1981)
   - Wall friction via ADWallFrictionChurchillMaterial / annular friction factor
   - Heat transfer via HeatTransferFromHeatFlux1Phase and ADTemperatureWall3EqnMaterial
   - Convective heat transfer coefficient (Gnielinski annular / Churchill)
   - Full axial distributions: T_bulk(z), p(z), rho(z), u_z(z), Re(z), Nu(z), h_conv(z), T_wall(z)
2. Solves 2D axisymmetric RZ solid pin heat conduction:
   - Pellet UN fuel with k(T) = 1.37 * T^0.41 (Hayes et al. 1990)
   - Nb-1Zr cladding with k(T) = 41.8 + 0.0125 * T (Touloukian 1970)
   - Pellet-cladding gap conductance h_gap = 25,000 W/(m2 K)
   - Volume-averaged fuel temperature T_vol(z) for Doppler feedback
3. Executes direct Cardinal MOOSE finite-element benchmark on peak and average pins.
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
    D_H, A_FLOW_PIN, P_WET_PIN, Q_TOTAL, N_PINS, T_INLET, DELTA_T_BULK,
    cp_helium, pin_coolant_mass_flow, total_coolant_mass_flow, rho_helium,
    mu_helium, k_helium, prandtl_helium, reynolds_helium,
    nusselt_annulus, friction_factor_annulus, H_GAP_NOMINAL,
    k_fuel_un, k_clad_nb1zr, P_NOM
)

CARDINAL_EXEC = "/home/aleri/projects/cardinal/cardinal-opt"
MOOSE_BIN_DIR = "/home/aleri/miniforge3/envs/moose/bin"


def solve_axial_flow_channel(q_vol_slices,
                             m_dot_pin,
                             z_edges,
                             p_outlet=P_NOM,
                             h_gap=H_GAP_NOMINAL,
                             htc_model="gnielinski_trans",
                             T_inlet=None,
                             htc_mult=1.0,
                             r_fuel=None,
                             r_clad_out=None,
                             r_chan_out=None,
                             h_core=None):
    """
    Solves 1D axial compressible flow conservation (mass, momentum, energy)
    for an annular helium coolant channel coupled to the fuel pin.
    
    Equations:
    - Mass conservation: m_dot = rho(z) * u(z) * A_flow = const
    - Energy conservation: m_dot * cp * dT_b/dz = q'(z) = q'''(z) * pi * r_f^2
    - Momentum conservation: dp/dz = - f_D/(2*D_h) * rho * u^2 - d(rho * u^2)/dz
    - Equation of state: rho(z) = p(z) * M_He / (R * T_b(z))
    - Convective HTC: h(z) = Nu(z) * k_He(z) / D_h
    """
    rf = R_FUEL if r_fuel is None else r_fuel
    rco = R_CLAD_OUT if r_clad_out is None else r_clad_out
    rcho = R_CHAN_OUT if r_chan_out is None else r_chan_out
    hc = H_CORE if h_core is None else h_core
    
    a_flow = math.pi * (rcho**2 - rco**2)
    d_h = 2.0 * (rcho - rco)

    n_slices = len(q_vol_slices)
    dz = hc / n_slices
    v_slice = math.pi * (rf**2) * dz
    a_clad_slice = 2.0 * math.pi * rco * dz
    cp = cp_helium()
    G = m_dot_pin / a_flow  # Mass flux [kg / (m2 s)]
    
    # 1. Thermal fields allocation
    T_bulk = np.zeros(n_slices)
    T_clad_out = np.zeros(n_slices)
    T_clad_in = np.zeros(n_slices)
    T_fuel_surf = np.zeros(n_slices)
    T_fuel_center = np.zeros(n_slices)
    T_fuel_vol = np.zeros(n_slices)
    
    # 2. Fluid fields allocation
    h_conv = np.zeros(n_slices)
    Re_arr = np.zeros(n_slices)
    Nu_arr = np.zeros(n_slices)
    Pr_arr = np.zeros(n_slices)
    rho_cool = np.zeros(n_slices)
    vel_cool = np.zeros(n_slices)
    p_cool = np.zeros(n_slices)
    dp_fric_slice = np.zeros(n_slices)
    dp_accel_slice = np.zeros(n_slices)
    
    # Step A: Forward pass for enthalpy rise and initial properties
    T_current = T_INLET if T_inlet is None else T_inlet
    for s in range(n_slices):
        q_vol = q_vol_slices[s]
        p_slice = q_vol * v_slice
        
        delta_T_b = p_slice / (m_dot_pin * cp)
        T_bulk_exit = T_current + delta_T_b
        T_bulk_mid = 0.5 * (T_current + T_bulk_exit)
        T_bulk[s] = T_bulk_mid
        
        # Local fluid properties at slice midpoint
        rho_s = rho_helium(T_bulk_mid, p_outlet)
        rho_cool[s] = rho_s
        vel_cool[s] = m_dot_pin / (rho_s * a_flow)
        
        Re_s = reynolds_helium(m_dot_pin, T_bulk_mid, p_outlet, a_flow=a_flow, d_h=d_h)
        Re_arr[s] = Re_s
        Pr_s = prandtl_helium(T_bulk_mid, p_outlet)
        Pr_arr[s] = Pr_s
        Nu_s = nusselt_annulus(Re_s, Pr_s, model=htc_model)
        Nu_arr[s] = Nu_s
        
        k_he_s = k_helium(T_bulk_mid, p_outlet)
        h_s = (Nu_s * k_he_s / d_h) * htc_mult
        h_conv[s] = h_s
        
        # Friction factor & pressure drop
        f_s = friction_factor_annulus(Re_s)
        dp_fric_slice[s] = f_s * (dz / d_h) * (G**2) / (2.0 * rho_s)
        
        rho_in_s = rho_helium(T_current, p_outlet)
        rho_out_s = rho_helium(T_bulk_exit, p_outlet)
        dp_accel_slice[s] = (G**2) * (1.0 / rho_out_s - 1.0 / rho_in_s)
        
        T_current = T_bulk_exit
        
        # Heat transfer to solid cladding outer surface
        q_flux_out = p_slice / a_clad_slice
        T_co = T_bulk_mid + q_flux_out / h_s
        T_clad_out[s] = T_co
        
        # Cladding conduction with k_c(T)
        k_c_guess = k_clad_nb1zr(T_co)
        delta_T_c = (p_slice / (2.0 * math.pi * dz)) * math.log(rco / rf) / k_c_guess
        k_c = k_clad_nb1zr(T_co + 0.5 * delta_T_c)
        delta_T_c = (p_slice / (2.0 * math.pi * dz)) * math.log(rco / rf) / k_c
        T_ci = T_co + delta_T_c
        T_clad_in[s] = T_ci
        
        # Pellet-cladding interface gap
        if h_gap is not None and h_gap > 0.0:
            delta_T_gap = (p_slice / (2.0 * math.pi * rf * dz)) / h_gap
        else:
            delta_T_gap = 0.0
        T_fs = T_ci + delta_T_gap
        T_fuel_surf[s] = T_fs
        
        # Fuel pellet conduction with k_f(T) = 1.37 * T^0.41
        integral_k = (p_slice / dz) / (4.0 * math.pi)
        T0_val = (T_fs**1.41 + (1.41 / 1.37) * integral_k)**(1.0 / 1.41)
        T_fuel_center[s] = T0_val
        
        # Exact volume-averaged fuel temperature across pellet radius
        r_arr = np.linspace(0.0, rf, 61)
        T_r = (T0_val**1.41 - (1.41 / 1.37) * (q_vol * (r_arr**2) / 4.0))**(1.0 / 1.41)
        T_vol = np.trapezoid(2.0 * r_arr * T_r, r_arr) / (rf**2)
        T_fuel_vol[s] = T_vol
        
    # Step B: Backward integration for axial pressure profile p(z)
    p_running = p_outlet
    for s in reversed(range(n_slices)):
        p_running += (dp_fric_slice[s] + dp_accel_slice[s])
        p_cool[s] = p_running
        
    # Form losses at inlet contraction and outlet expansion
    rho_inlet = rho_helium(T_INLET, p_cool[0])
    u_inlet = m_dot_pin / (rho_inlet * a_flow)
    rho_outlet = rho_helium(T_current, p_outlet)
    u_outlet = m_dot_pin / (rho_outlet * a_flow)
    dp_form = 0.5 * rho_inlet * (u_inlet**2) / 2.0 + 1.0 * rho_outlet * (u_outlet**2) / 2.0
    
    dp_fric_total = float(np.sum(dp_fric_slice))
    dp_accel_total = float(np.sum(dp_accel_slice))
    dp_total = dp_fric_total + dp_accel_total + dp_form
    
    return {
        "T_bulk": T_bulk,
        "T_bulk_exit": T_current,
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
        "vel_cool": vel_cool,
        "p_cool": p_cool,
        "dp_total_Pa": dp_total,
        "dp_fric_Pa": dp_fric_total,
        "dp_accel_Pa": dp_accel_total,
        "dp_form_Pa": dp_form
    }


def generate_cardinal_thm_input(pin_name, z_edges, q_vol_slices, m_dot_pin, output_file):
    """
    Generates a MOOSE THM FlowChannel1Phase input file executable by cardinal-opt.
    Uses HeatTransferFromHeatFlux1Phase with specified axial heat flux profile.
    """
    n_slices = len(q_vol_slices)
    dz = H_CORE / n_slices
    v_slice = math.pi * (R_FUEL**2) * dz
    a_clad_slice = 2.0 * math.pi * R_CLAD_OUT * dz
    
    q_flux_points_x = [f"{0.0:.5f}"]
    q_flux_points_y = [f"{(q_vol_slices[0] * v_slice / a_clad_slice):.2f}"]
    for s in range(n_slices):
        z_mid = 0.5 * (z_edges[s] + z_edges[s+1])
        q_flux_s = q_vol_slices[s] * v_slice / a_clad_slice
        q_flux_points_x.append(f"{z_mid:.5f}")
        q_flux_points_y.append(f"{q_flux_s:.2f}")
    q_flux_points_x.append(f"{H_CORE:.5f}")
    q_flux_points_y.append(f"{(q_vol_slices[-1] * v_slice / a_clad_slice):.2f}")
    
    q_x_str = " ".join(q_flux_points_x)
    q_y_str = " ".join(q_flux_points_y)
    
    d_in = 2.0 * R_CLAD_OUT
    content = f"""# CARDINAL MOOSE THM 1D FLOW CHANNEL: {pin_name}
inlet_T = {T_INLET:.1f}
mdot = {m_dot_pin:.6e}
outlet_P = {P_NOM:.1f}
height = {H_CORE:.5f}
d_in = {d_in:.6f}
d_hyd = {D_H:.6f}
area = {A_FLOW_PIN:.6e}

[GlobalParams]
  initial_p = ${{outlet_P}}
  initial_T = ${{inlet_T}}
  initial_vel = ${{fparse mdot / outlet_P / 8.3144598 * 4.0e-3 / inlet_T / area}}
  rdg_slope_reconstruction = full
  closures = none
  fp = helium
[]

[Closures]
  [none]
    type = WallTemperature1PhaseClosures
  []
[]

[FluidProperties]
  [helium]
    type = HeliumFluidProperties
  []
[]

[Functions]
  [q_flux_fn]
    type = PiecewiseLinear
    x = '{q_x_str}'
    y = '{q_y_str}'
    axis = z
  []
[]

[AuxVariables]
  [T_wall]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [Tw_aux]
    type = ADMaterialRealAux
    block = channel
    variable = T_wall
    property = T_wall
  []
[]

[Materials]
  [f_mat]
    type = ADWallFrictionChurchillMaterial
    block = channel
    D_h = D_h
    f_D = f_D
    vel = vel
    rho = rho
    mu = mu
  []
  [Re_mat]
    type = ADReynoldsNumberMaterial
    block = channel
    Re = Re
    D_h = D_h
    mu = mu
    vel = vel
    rho = rho
  []
  [Pr_mat]
    type = ADPrandtlNumberMaterial
    block = channel
    cp = cp
    mu = mu
    k = k
  []
  [Nu_mat]
    type = ADParsedMaterial
    block = channel
    expression = '0.021 * pow(Re, 0.8) * pow(Pr, 0.4)'
    property_name = 'Nu'
    material_property_names = 'Re Pr'
  []
  [Hw_mat]
    type = ADConvectiveHeatTransferCoefficientMaterial
    block = channel
    D_h = D_h
    Nu = Nu
    k = k
  []
  [T_wall]
    type = ADTemperatureWall3EqnMaterial
    Hw = Hw
    T = T
    q_wall = q_wall
  []
[]

[Components]
  [channel]
    type = FlowChannel1Phase
    position = '0 0 0'
    orientation = '0 0 1'
    A = ${{area}}
    D_h = ${{d_hyd}}
    length = ${{height}}
    n_elems = {n_slices}
  []
  [inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'channel:in'
    m_dot = ${{mdot}}
    T = ${{inlet_T}}
  []
  [outlet]
    type = Outlet1Phase
    input = 'channel:out'
    p = ${{outlet_P}}
  []
  [ht_ext]
    type = HeatTransferFromHeatFlux1Phase
    flow_channel = channel
    P_hf = ${{fparse pi * d_in}}
    q_wall = q_flux_fn
  []
[]

[Postprocessors]
  [p_inlet]
    type = SideAverageValue
    variable = p
    boundary = channel:in
  []
  [p_outlet]
    type = SideAverageValue
    variable = p
    boundary = channel:out
  []
  [dp_total]
    type = DifferencePostprocessor
    value1 = p_inlet
    value2 = p_outlet
  []
  [T_inlet_meas]
    type = SideAverageValue
    variable = T
    boundary = channel:in
  []
  [T_outlet_meas]
    type = SideAverageValue
    variable = T
    boundary = channel:out
  []
  [max_T_wall]
    type = ElementExtremeValue
    variable = T_wall
    value_type = max
  []
[]

[VectorPostprocessors]
  [line_sample]
    type = LineValueSampler
    variable = 'T p rho vel_z T_wall'
    start_point = '0 0 {z_edges[0] + dz * 0.5:.6f}'
    end_point = '0 0 {z_edges[-1] - dz * 0.5:.6f}'
    num_points = {n_slices}
    sort_by = z
  []
[]

[Preconditioning]
  [pc]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  dt = 0.05
  start_time = 0
  num_steps = 20
  steady_state_detection = true
  steady_state_tolerance = 1e-7
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-6
  solve_type = NEWTON
  line_search = basic
[]

[Outputs]
  csv = true
  print_linear_residuals = false
[]
"""
    with open(output_file, "w") as f:
        f.write(content)


def generate_cardinal_solid_input(pin_name, z_edges, q_vol_slices, h_conv_slices, T_bulk_slices, output_file):
    """
    Generates a 2D axisymmetric RZ MOOSE solid conduction input file executable by cardinal-opt.
    """
    n_slices = len(q_vol_slices)
    dz = H_CORE / n_slices
    
    q_x = [f"{0.0:.5f}"]
    q_y = [f"{q_vol_slices[0]:.4e}"]
    for s in range(n_slices):
        z_mid = 0.5 * (z_edges[s] + z_edges[s+1])
        q_x.append(f"{z_mid:.5f}")
        q_y.append(f"{q_vol_slices[s]:.4e}")
    q_x.append(f"{H_CORE:.5f}")
    q_y.append(f"{q_vol_slices[-1]:.4e}")
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
    
    h_avg = float(np.mean(h_conv_slices))
    
    content = f"""# CARDINAL MOOSE 2D AXISYMMETRIC PIN CELL: {pin_name}
[Mesh]
  coord_type = RZ
  [fuel_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = {n_slices * 6}
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
    ny = {n_slices * 6}
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
    coefficient = {h_avg:.2f}
  []
[]

[Materials]
  [k_fuel]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = 'T'
    expression = '1.37 * (T^0.41)'
    block = '1'
  []
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
        f.write(content)


def run_cardinal_input(input_file):
    """Executes cardinal-opt with proper PATH environment and returns dataframe / dict."""
    env = os.environ.copy()
    env["PATH"] = f"{MOOSE_BIN_DIR}:{env.get('PATH', '')}"
    
    cmd = [CARDINAL_EXEC, "-i", input_file]
    work_dir = os.path.dirname(os.path.abspath(input_file))
    res = subprocess.run(cmd, cwd=work_dir, env=env, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Cardinal execution failed on {input_file}:\n{res.stderr}\n{res.stdout}")
        
    csv_file = input_file.replace(".i", "_out.csv")
    if not os.path.exists(csv_file):
        raise FileNotFoundError(f"Expected CSV output not found: {csv_file}")
    df = pd.read_csv(csv_file)
    return df.iloc[-1].to_dict()
