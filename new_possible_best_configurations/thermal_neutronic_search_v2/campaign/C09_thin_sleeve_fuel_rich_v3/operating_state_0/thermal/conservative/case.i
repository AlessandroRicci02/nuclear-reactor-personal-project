[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010795
    xmax = 0.012295
    ymin = 0
    ymax = 0.72
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
    x = '0 0.09 0.18 0.27 0.36 0.45 0.54 0.63 0.72'
    y = '2295270.30657 2434141.27836 3090709.48833 3437102.93734 3442764.51873 3105438.90687 2450563.52977 2301296.66914 2301296.66914'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.09 0.18 0.27 0.36 0.45 0.54 0.63 0.72'
    y = '1005.39984717 993.844941158 983.208926139 968.905439971 951.418661144 931.523556905 910.280283257 889.006753525 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.09 0.18 0.27 0.36 0.45 0.54 0.63 0.72'
    y = '1005.39984717 993.844941158 983.208926139 968.905439971 951.418661144 931.523556905 910.280283257 889.006753525 850'
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
    coefficient = 564.744656638578
    T_infinity = inner_bulk
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
