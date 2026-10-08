[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010395000000000001
    xmax = 0.012395
    ymin = 0
    ymax = 0.7
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
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '1000483.19357 1087754.50907 1394956.16496 1545205.58989 1550866.07908 1392869.82233 1100536.74464 1003401.66834 1003401.66834'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '930.10675553 927.777018209 925.328887325 922.108024614 918.372782615 914.385310394 910.529959836 907.643694662 900'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '930.227278199 927.909508912 925.332835838 921.854937219 917.761827104 913.324560582 908.820005805 903.942647344 900'
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
    coefficient = 798.6230744212256
    T_infinity = inner_bulk
  []
  [outer]
    type = ConvectiveFluxFunction
    variable = T
    boundary = right
    coefficient = 814.2478978143457
    T_infinity = outer_bulk
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
