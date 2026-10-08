[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010395000000000001
    xmax = 0.012895
    ymin = 0
    ymax = 0.68
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
    x = '0 0.085 0.17 0.255 0.34 0.425 0.51 0.595 0.68'
    y = '1518324.8321 1475958.23379 1842188.91492 2042312.12635 2040888.82038 1843822.92105 1465566.97111 1510841.54784 1510841.54784'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.085 0.17 0.255 0.34 0.425 0.51 0.595 0.68'
    y = '1019.71038395 1004.65517171 992.026240078 975.777189977 956.552538043 935.436324737 912.849517253 895.107145877 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.085 0.17 0.255 0.34 0.425 0.51 0.595 0.68'
    y = '1020.03904376 1006.02942733 993.444380842 976.659112091 956.348687765 933.398356629 910.683200831 879.588850778 850'
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
    coefficient = 645.0780629439561
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
