[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010495
    xmax = 0.012495000000000001
    ymin = 0
    ymax = 0.8
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
    x = '0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8'
    y = '1158794.64883 1278667.38965 1639365.49419 1826324.96195 1828517.88629 1630177.68324 1275109.15059 1146919.91065 1146919.91065'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8'
    y = '921.148497602 917.733726562 913.942389725 908.867694295 902.873847605 896.389726553 889.93810601 883.982036454 875'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8'
    y = '921.148497602 917.733726562 913.942389725 908.867694295 902.873847605 896.389726553 889.93810601 883.982036454 875'
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
