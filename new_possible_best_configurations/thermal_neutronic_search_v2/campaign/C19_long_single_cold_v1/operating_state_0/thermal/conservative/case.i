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
    y = '2310179.3182 2470159.61019 3118959.60324 3447836.40397 3440389.35895 3080399.31204 2434538.71469 2280928.23694 2280928.23694'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.095 0.19 0.285 0.38 0.475 0.57 0.665 0.76'
    y = '978.613030549 968.6599966 959.002739652 946.234597979 930.970954343 914.000481255 896.371797464 879.180282232 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.095 0.19 0.285 0.38 0.475 0.57 0.665 0.76'
    y = '978.613030549 968.6599966 959.002739652 946.234597979 930.970954343 914.000481255 896.371797464 879.180282232 850'
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
    coefficient = 564.3970755865249
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
