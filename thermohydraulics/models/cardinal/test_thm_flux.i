inlet_T = 950.0
mdot = 1.774691e-03
outlet_P = 2.5e6
height = 0.4245
d_in = 0.0126
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

[Functions]
  [q_flux_fn]
    type = PiecewiseLinear
    x = '0.0 0.04245 0.08490 0.12735 0.16980 0.21225 0.25470 0.29715 0.33960 0.38205 0.42450'
    y = '34000 50000 70000 90000 100000 102000 95000 80000 60000 40000 30000'
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
    n_elems = 10
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
    type = HeatTransferFromHeatFlux1Phase
    flow_channel = channel
    P_hf = ${fparse pi * d_in}
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
    start_point = '0 0 0.021225'
    end_point = '0 0 0.403275'
    num_points = 10
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
