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
    y = '561847.889576 561847.889576 561847.889576 561847.889576 561847.889576 561847.889576 561847.889576 561847.889576 561847.889576'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0775 0.155 0.2325 0.31 0.3875 0.465 0.5425 0.62'
    y = '875.826924549 874.039509569 872.328080758 870.333784698 868.030730361 865.400134178 862.518212598 860.569400932 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0775 0.155 0.2325 0.31 0.3875 0.465 0.5425 0.62'
    y = '875.629514874 873.765476242 871.717307509 869.355382087 866.650873318 863.564529493 859.947909065 854.401492084 850'
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
