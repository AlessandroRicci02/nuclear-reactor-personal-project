[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010795
    xmax = 0.011995
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
    y = '799178.673166 799178.673166 799178.673166 799178.673166 799178.673166 799178.673166 799178.673166 799178.673166 799178.673166'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0825 0.165 0.2475 0.33 0.4125 0.495 0.5775 0.66'
    y = '880.279586168 878.326601789 876.308120189 873.894773095 871.049914003 867.73583344 863.914058398 859.545687103 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0825 0.165 0.2475 0.33 0.4125 0.495 0.5775 0.66'
    y = '880.279586168 878.326601789 876.308120189 873.894773095 871.049914003 867.73583344 863.914058398 859.545687103 850'
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
    coefficient = 799.2822414908338
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
