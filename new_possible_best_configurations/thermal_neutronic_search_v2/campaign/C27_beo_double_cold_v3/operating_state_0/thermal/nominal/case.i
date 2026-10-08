[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010595
    xmax = 0.012895
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
    y = '733151.125991 803282.025447 1022054.93106 1143666.96791 1131054.68431 1029157.1261 802323.429307 731556.032678 731556.032678'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '917.050539372 912.600102666 907.811890181 901.367818506 893.52937844 885.099529791 875.343507168 869.749449546 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '917.026740842 912.60584558 907.606561577 900.775599613 892.509177966 883.110705275 874.332162922 861.516148378 850'
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
