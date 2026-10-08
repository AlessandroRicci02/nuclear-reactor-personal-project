[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010795
    xmax = 0.012595
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
    y = '2590091.23256 2685737.54839 3428092.58252 3787795.9894 3800420.32048 3397154.99653 2702032.75395 2594505.83621 2594505.83621'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0825 0.165 0.2475 0.33 0.4125 0.495 0.5775 0.66'
    y = '953.950127469 945.213451952 937.159172067 926.561169746 914.086336234 900.42450806 886.533857296 873.226635337 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0825 0.165 0.2475 0.33 0.4125 0.495 0.5775 0.66'
    y = '953.950127469 945.213451952 937.159172067 926.561169746 914.086336234 900.42450806 886.533857296 873.226635337 850'
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
