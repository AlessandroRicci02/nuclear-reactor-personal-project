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
    y = '596092.560804 584433.8347 740505.65125 828452.954796 825910.383123 743413.349318 577943.243653 588017.719689 588017.719689'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.085 0.17 0.255 0.34 0.425 0.51 0.595 0.68'
    y = '919.804388407 915.018270092 910.317279731 903.977373859 896.074356584 887.6896823 877.294062562 871.76998237 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.085 0.17 0.255 0.34 0.425 0.51 0.595 0.68'
    y = '919.734878382 914.947343458 910.039452112 903.265340165 895.122484863 885.40740417 876.846217942 863.235332213 850'
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
    coefficient = 813.5872529931196
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
