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
    y = '2091299.10803 2178551.91934 2747497.67241 3007895.93938 2995166.18991 2700556.71934 2121366.82235 2041426.82333 2041426.82333'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.105 0.21 0.315 0.42 0.525 0.63 0.735 0.84'
    y = '992.044589303 983.992305337 976.539541882 966.890187261 955.64795093 943.423778054 930.957847675 919.112661584 900'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.105 0.21 0.315 0.42 0.525 0.63 0.735 0.84'
    y = '992.044589303 983.992305337 976.539541882 966.890187261 955.64795093 943.423778054 930.957847675 919.112661584 900'
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
    coefficient = 564.3970755865249
    T_infinity = inner_bulk
  []
[]
[Materials/thermal]
  type = GenericConstantMaterial
  prop_names = thermal_conductivity
  prop_values = 5.0
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
