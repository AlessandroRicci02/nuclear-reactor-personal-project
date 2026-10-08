[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010395000000000001
    xmax = 0.012595
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
    y = '10646695.0066 7241796.8602 8708302.20941 8605970.22062 10516605.4147 12493260.6957 8112936.31473 9835185.6772 9835185.6772'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0775 0.155 0.2325 0.31 0.3875 0.465 0.5425 0.62'
    y = '1189.84843218 1152.81113457 1122.3519404 1085.35418389 1045.07817977 996.921545541 943.443791283 899.04858641 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0775 0.155 0.2325 0.31 0.3875 0.465 0.5425 0.62'
    y = '1209.80637782 1164.22828439 1134.90254745 1096.13489652 1058.64837186 1012.36483619 948.415274438 906.099601182 850'
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
    coefficient = 564.2767561897358
    T_infinity = inner_bulk
  []
  [outer]
    type = ConvectiveFluxFunction
    variable = T
    boundary = right
    coefficient = 645.4144161358766
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
