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
    ymax = 0.8
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
    x = '0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8'
    y = '2105106.97615 2221310.6065 2810760.16415 3100713.07949 3116379.82344 2789662.42075 2211581.60361 2087276.48298 2087276.48298'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8'
    y = '975.7061432 967.470310962 959.635911522 949.373207304 937.236836058 923.827950362 910.06942106 896.788738131 875'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8'
    y = '975.7061432 967.470310962 959.635911522 949.373207304 937.236836058 923.827950362 910.06942106 896.788738131 875'
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
