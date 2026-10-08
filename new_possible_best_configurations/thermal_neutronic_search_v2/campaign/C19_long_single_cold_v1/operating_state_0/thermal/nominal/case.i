[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010495
    xmax = 0.012495000000000001
    ymin = 0
    ymax = 0.76
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
    x = '0 0.095 0.19 0.285 0.38 0.475 0.57 0.665 0.76'
    y = '959396.105417 1083529.73391 1380696.5025 1542044.56705 1537564.36283 1371318.69868 1068944.88814 947016.65765 947016.65765'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.095 0.19 0.285 0.38 0.475 0.57 0.665 0.76'
    y = '897.83048137 894.645941509 891.014888804 886.109019712 880.199851551 873.676736402 866.985312594 860.584984505 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.095 0.19 0.285 0.38 0.475 0.57 0.665 0.76'
    y = '897.83048137 894.645941509 891.014888804 886.109019712 880.199851551 873.676736402 866.985312594 860.584984505 850'
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
    coefficient = 798.7925778545774
    T_infinity = inner_bulk
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
