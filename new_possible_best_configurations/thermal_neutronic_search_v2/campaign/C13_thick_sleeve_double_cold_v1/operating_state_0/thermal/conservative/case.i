[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010195000000000001
    xmax = 0.013195000000000002
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
    y = '2264985.56799 2334819.71861 2943987.96263 3253372.82985 3242784.60828 2932348.67667 2320632.06336 2238613.60056 2238613.60056'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0775 0.155 0.2325 0.31 0.3875 0.465 0.5425 0.62'
    y = '964.529098803 954.271528981 944.742662282 932.768785296 919.096645151 904.625666026 890.444834262 878.217656747 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0775 0.155 0.2325 0.31 0.3875 0.465 0.5425 0.62'
    y = '968.412035096 959.122304656 949.97626796 937.250823715 921.759429074 904.491732298 886.243055734 867.054822429 850'
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
    coefficient = 564.0290351222121
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
