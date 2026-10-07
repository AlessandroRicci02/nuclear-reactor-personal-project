# MOOSE 2D Axisymmetric FEM Pin Cell Model: Peak_Pin_203
[Mesh]
  coord_type = RZ
  [fuel_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 60
    xmin = 0.0
    xmax = 0.00580
    ymin = 0.0
    ymax = 0.42450
    boundary_name_prefix = fuel
  []
  [fuel_rename]
    type = SubdomainIDGenerator
    input = fuel_mesh
    subdomain_id = 1
  []
  [clad_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 60
    xmin = 0.00580
    xmax = 0.00628
    ymin = 0.0
    ymax = 0.42450
    boundary_name_prefix = clad
  []
  [clad_rename]
    type = SubdomainIDGenerator
    input = clad_mesh
    subdomain_id = 2
  []
  [stitch]
    type = StitchMeshGenerator
    inputs = 'fuel_rename clad_rename'
    stitch_boundaries_pairs = 'fuel_right clad_left'
  []
[]

[Variables]
  [T]
    initial_condition = 950.0
  []
[]

[Functions]
  [q3_func]
    type = PiecewiseLinear
    x = '0.00000 0.02123 0.06368 0.10613 0.14858 0.19103 0.23348 0.27593 0.31838 0.36083 0.40327 0.42450'
    y = '3.0254e+07 3.0254e+07 4.5474e+07 5.5196e+07 6.3888e+07 6.8405e+07 6.9386e+07 6.8433e+07 5.5885e+07 4.5011e+07 3.0624e+07 3.0624e+07'
    axis = y
  []
  [T_bulk_func]
    type = PiecewiseLinear
    x = '0.00000 0.02123 0.06368 0.10613 0.14858 0.19103 0.23348 0.27593 0.31838 0.36083 0.40327 0.42450'
    y = '950.000 961.486 978.749 999.704 1023.958 1049.927 1076.269 1102.249 1123.465 1140.553 1152.179 1152.179'
    axis = y
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
  []
  [source]
    type = BodyForce
    variable = T
    function = q3_func
    block = '1'
  []
[]

[BCs]
  [outer_convection]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'clad_right'
    T_infinity = T_bulk_func
    coefficient = 1229.71
  []
[]

[Materials]
  # UN Fuel: k = 1.37 * T^0.41 [Hayes & Peddicord 1990]
  [k_fuel]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = 'T'
    expression = '1.37 * (T^0.41)'
    block = '1'
  []
  # Nb-1Zr Clad: k = 42.0 + 0.015 * (T - 293.15) [TPRC Series]
  [k_clad]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = 'T'
    expression = '42.0 + 0.015 * (T - 293.15)'
    block = '2'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-6
  nl_max_its = 30
  l_max_its = 100
[]

[Postprocessors]
  [max_T_fuel]
    type = ElementExtremeValue
    variable = T
    value_type = max
    block = '1'
  []
  [avg_T_fuel]
    type = ElementAverageValue
    variable = T
    block = '1'
  []
  [max_T_clad]
    type = ElementExtremeValue
    variable = T
    value_type = max
    block = '2'
  []
  [avg_T_clad]
    type = ElementAverageValue
    variable = T
    block = '2'
  []
  [max_T_clad_surf]
    type = SideExtremeValue
    variable = T
    boundary = 'clad_right'
    value_type = max
  []
  [avg_T_clad_surf]
    type = SideAverageValue
    variable = T
    boundary = 'clad_right'
  []
[]

[Outputs]
  csv = true
  console = false
[]
