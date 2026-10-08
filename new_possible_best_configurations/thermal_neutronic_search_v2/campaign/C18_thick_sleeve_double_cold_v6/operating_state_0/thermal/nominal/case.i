[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010195000000000001
    xmax = 0.012995000000000001
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
    y = '640289.571773 651842.691634 828050.524305 924231.835295 926350.410588 830117.601332 651393.445773 640325.034442 640325.034442'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '935.235742119 932.409837208 929.601323224 925.913592738 921.583368106 916.931010191 912.326192793 909.186095104 900'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '935.301217859 932.477807225 929.53460954 925.600955761 920.935880801 915.887331144 910.849032892 905.021760732 900'
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
    coefficient = 798.2740895046661
    T_infinity = inner_bulk
  []
  [outer]
    type = ConvectiveFluxFunction
    variable = T
    boundary = right
    coefficient = 813.4611876553071
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
