[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010495
    xmax = 0.013095
    ymin = 0
    ymax = 0.76
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
    x = '0 0.095 0.19 0.285 0.38 0.475 0.57 0.665 0.76'
    y = '506986.694303 514277.573278 655671.204128 737136.011836 732198.417902 658807.919873 510912.772202 510402.274913 510402.274913'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.095 0.19 0.285 0.38 0.475 0.57 0.665 0.76'
    y = '898.845385162 895.364114492 891.776821738 886.912759833 881.00793014 874.481502318 867.716521184 861.217916093 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.095 0.19 0.285 0.38 0.475 0.57 0.665 0.76'
    y = '898.845385162 895.364114492 891.776821738 886.912759833 881.00793014 874.481502318 867.716521184 861.217916093 850'
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
