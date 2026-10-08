[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010195000000000001
    xmax = 0.013795000000000002
    ymin = 0
    ymax = 0.62
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
    x = '0 0.0775 0.155 0.2325 0.31 0.3875 0.465 0.5425 0.62'
    y = '1540499.4406 1475812.39334 1841839.61396 2041921.59259 2042627.73272 1840056.49838 1472783.35718 1564030.93034 1564030.93034'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0775 0.155 0.2325 0.31 0.3875 0.465 0.5425 0.62'
    y = '947.856517323 938.042545513 929.541975627 919.055679647 907.203168497 894.776546554 882.751887353 872.512692806 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0775 0.155 0.2325 0.31 0.3875 0.465 0.5425 0.62'
    y = '950.724954774 941.65855867 933.602887317 922.774573469 909.759171897 895.41217252 880.598582714 865.18731989 850'
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
    coefficient = 564.0290351222121
    T_infinity = inner_bulk
  []
  [outer]
    type = ConvectiveFluxFunction
    variable = T
    boundary = right
    coefficient = 644.1563693241131
    T_infinity = outer_bulk
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
