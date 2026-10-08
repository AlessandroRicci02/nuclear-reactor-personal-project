[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010495
    xmax = 0.012295
    ymin = 0
    ymax = 0.84
  []
[]
[Variables/T]
  initial_condition = 900
[]
[Functions]
  [source]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '0 0.105 0.21 0.315 0.42 0.525 0.63 0.735 0.84'
    y = '861893.617488 950011.51265 1205659.51039 1329929.8795 1320703.27345 1185697.60744 915023.081175 841414.545597 841414.545597'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.105 0.21 0.315 0.42 0.525 0.63 0.735 0.84'
    y = '934.584243668 932.045499234 929.243551171 925.503295944 921.091354453 916.310018542 911.480311032 906.994446011 900'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.105 0.21 0.315 0.42 0.525 0.63 0.735 0.84'
    y = '934.584243668 932.045499234 929.243551171 925.503295944 921.091354453 916.310018542 911.480311032 906.994446011 900'
  []
[]
[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
  []
  [heat]
    type = BodyForce
    variable = T
    function = source
  []
[]
[BCs]
  [inner]
    type = ConvectiveFluxFunction
    variable = T
    boundary = left
    coefficient = 798.7925778545774
    T_infinity = inner_bulk
  []
[]
[Materials/thermal]
  type = GenericConstantMaterial
  prop_names = thermal_conductivity
  prop_values = 10.0
[]
[Postprocessors/Tmax]
  type = NodalExtremeValue
  variable = T
  value_type = max
[]
[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
[]
[Outputs]
  csv = true
  exodus = true
[]
