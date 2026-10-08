[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010795
    xmax = 0.012295
    ymin = 0
    ymax = 0.72
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
    x = '0 0.09 0.18 0.27 0.36 0.45 0.54 0.63 0.72'
    y = '1052054.65341 1158622.04149 1471196.20004 1650296.5737 1650089.67957 1487306.17708 1154896.80567 1057548.79153 1057548.79153'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.09 0.18 0.27 0.36 0.45 0.54 0.63 0.72'
    y = '905.08781767 901.614369647 897.755842247 892.48617769 886.026171145 878.765596657 871.111518682 863.602635784 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.09 0.18 0.27 0.36 0.45 0.54 0.63 0.72'
    y = '905.08781767 901.614369647 897.755842247 892.48617769 886.026171145 878.765596657 871.111518682 863.602635784 850'
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
    coefficient = 799.2822414908338
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
