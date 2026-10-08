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
    y = '1518860.25863 1102764.93978 1026353.32465 1250172.45975 1549213.93835 1951522.46977 1281143.63681 1287767.43653 1287767.43653'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0775 0.155 0.2325 0.31 0.3875 0.465 0.5425 0.62'
    y = '906.656296695 900.905137142 896.451035935 891.62077175 885.44698462 877.773426171 868.094676351 860.30929943 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0775 0.155 0.2325 0.31 0.3875 0.465 0.5425 0.62'
    y = '907.326047256 901.147011403 896.219846738 891.156496318 884.831395788 876.872751356 866.34847212 857.597926266 850'
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
    coefficient = 813.9773848824279
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
