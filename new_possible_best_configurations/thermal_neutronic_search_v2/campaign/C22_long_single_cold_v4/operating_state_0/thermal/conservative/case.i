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
    y = '1321370.31355 1292878.38136 1623979.78263 1807334.99421 1801233.55801 1626005.40802 1281452.68693 1326958.4457 1326958.4457'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.095 0.19 0.285 0.38 0.475 0.57 0.665 0.76'
    y = '962.92628925 952.741028213 943.826733653 932.330935647 918.784835125 904.009650786 888.873647677 874.404056987 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.095 0.19 0.285 0.38 0.475 0.57 0.665 0.76'
    y = '962.92628925 952.741028213 943.826733653 932.330935647 918.784835125 904.009650786 888.873647677 874.404056987 850'
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
