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
    y = '1113599.14644 1189303.80601 1503081.40091 1668574.832 1662727.66205 1499800.13686 1182368.75948 1111130.0801 1111130.0801'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0825 0.165 0.2475 0.33 0.4125 0.495 0.5775 0.66'
    y = '920.193347299 916.473080267 912.558018755 907.511213267 901.72930047 895.684997506 889.929361446 885.601190391 875'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0825 0.165 0.2475 0.33 0.4125 0.495 0.5775 0.66'
    y = '920.548296845 916.812737137 912.741567499 907.362202172 901.076885173 894.381216526 887.669715298 880.785064069 875'
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
    coefficient = 798.2740895046661
    T_infinity = inner_bulk
  []
  [outer]
    type = ConvectiveFluxFunction
    variable = T
    boundary = right
    coefficient = 813.2147551541681
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
