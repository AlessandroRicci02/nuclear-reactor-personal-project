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
    y = '876341.912236 967415.747613 1240565.89482 1385046.90711 1383445.20345 1233551.12573 982874.772452 867172.437228 867172.437228'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.08 0.16 0.24 0.32 0.4 0.48 0.56 0.64'
    y = '893.558117382 890.256567533 886.705891554 882.026844746 876.564496147 870.735776795 865.045131338 861.108439388 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.08 0.16 0.24 0.32 0.4 0.48 0.56 0.64'
    y = '893.718195361 890.442857986 886.75546035 881.784807679 875.893926335 869.545422009 863.222578433 855.938085704 850'
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
    coefficient = 798.9588808694632
    T_infinity = inner_bulk
  []
  [outer]
    type = ConvectiveFluxFunction
    variable = T
    boundary = right
    coefficient = 813.2147551541681
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
