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
    ymax = 0.66
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
    x = '0 0.0825 0.165 0.2475 0.33 0.4125 0.495 0.5775 0.66'
    y = '2156678.59988 2227928.93599 2793587.50158 3088575.55127 3068998.74512 2785276.42162 2210132.55178 2152725.32096 2152725.32096'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0825 0.165 0.2475 0.33 0.4125 0.495 0.5775 0.66'
    y = '969.267137433 960.390793669 952.158831294 941.960927069 930.467888753 918.470421776 906.80054099 896.618460574 875'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0825 0.165 0.2475 0.33 0.4125 0.495 0.5775 0.66'
    y = '973.099236643 965.157616126 957.338169766 946.558977726 933.486650354 919.067496695 903.995984261 888.84369924 875'
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
