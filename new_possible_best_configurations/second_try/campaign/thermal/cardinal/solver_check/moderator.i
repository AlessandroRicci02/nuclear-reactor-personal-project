[Mesh]
  coord_type = RZ
  [moderator]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 12
    ny = 60
    xmin = 0.00927
    xmax = 0.01087
    ymin = 0
    ymax = 0.48
  []
[]
[Variables/T]
  initial_condition = 850
[]
[Functions]
  [heat]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '0 0.08 0.16 0.24 0.32 0.4 0.48'
    y = '1234754.709936 1234754.709936 1234754.709936 1234754.709936 1234754.709936 1234754.709936 1234754.709936'
  []
  [bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.08 0.16 0.24 0.32 0.4 0.48'
    y = '850 841.6666666667 833.3333333333 825 816.6666666667 808.3333333333 800'
  []
[]
[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
  []
  [source]
    type = BodyForce
    variable = T
    function = heat
  []
[]
[BCs]
  [inner]
    type = ConvectiveFluxFunction
    variable = T
    boundary = left
    coefficient = 1000.0
    T_infinity = bulk
  []
[]
[Materials]
  [conductivity]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 10.0
  []
  [bulk_temperature]
    type = GenericFunctionMaterial
    prop_names = bulk_temperature
    prop_values = bulk
  []
[]
[Postprocessors]
  [T_YH2_max_K]
    type = NodalExtremeValue
    variable = T
    value_type = max
  []
  [Q_removed_W]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    T_fluid = bulk_temperature
    htc = 1000.0
    boundary = left
  []
[]
[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-9
[]
[Outputs]
  csv = true
[]
