[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010595
    xmax = 0.013195000000000002
    ymin = 0
    ymax = 0.64
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
    x = '0 0.08 0.16 0.24 0.32 0.4 0.48 0.56 0.64'
    y = '2167393.05796 2261632.23876 2861909.59622 3184391.24576 3193152.50717 2853269.74512 2275360.56535 2158240.72678 2158240.72678'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.08 0.16 0.24 0.32 0.4 0.48 0.56 0.64'
    y = '958.674085323 948.659914887 939.427548508 927.814016145 914.526672432 900.467109246 886.816257439 875.311049224 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.08 0.16 0.24 0.32 0.4 0.48 0.56 0.64'
    y = '961.59674308 952.670512579 944.003691793 932.055697893 917.534342898 901.359316758 884.605540654 866.442489502 850'
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
    coefficient = 564.5151233670873
    T_infinity = inner_bulk
  []
  [outer]
    type = ConvectiveFluxFunction
    variable = T
    boundary = right
    coefficient = 644.7569296453586
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
