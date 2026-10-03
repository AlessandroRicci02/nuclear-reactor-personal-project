# Standalone MOOSE THM Flow Channel Test for Opt-17 Fuel Pin Annulus
inlet_T = 950.0
mdot = 1.774691e-03
outlet_P = 2.5e6
height = 0.4245
d_in = 0.0126
d_out = 0.0141
d_hyd = 0.0015
area = 3.14552e-05

[GlobalParams]
  initial_p = ${outlet_P}
  initial_T = ${inlet_T}
  initial_vel = ${fparse mdot / outlet_P / 8.3144598 * 4.0e-3 / inlet_T / area}
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
    A = ${area}
    D_h = ${d_hyd}
    length = ${height}
    n_elems = 20
  []
  [inlet]
    type = InletMassFlowRateTemperature1Phase
    input = 'channel:in'
    m_dot = ${mdot}
    T = ${inlet_T}
  []
  [outlet]
    type = Outlet1Phase
    input = 'channel:out'
    p = ${outlet_P}
  []
  [ht_ext]
    type = HeatTransferFromExternalAppHeatFlux1Phase
    flow_channel = channel
    P_hf = ${fparse pi * d_in}
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
  num_steps = 10
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
