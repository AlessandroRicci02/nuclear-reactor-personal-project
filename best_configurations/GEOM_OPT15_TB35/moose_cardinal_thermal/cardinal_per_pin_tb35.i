# MOOSE 2D Axisymmetric FEM Pin Cell Model: PERIPHERY_PIN_RING_8 (GEOM_OPT15_TB35)
[Mesh]
  coord_type = RZ
  [fuel_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 60
    xmin = 0.0
    xmax = 0.00585
    ymin = 0.0
    ymax = 0.45500
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
    xmin = 0.00585
    xmax = 0.00633
    ymin = 0.0
    ymax = 0.45500
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
    initial_condition = 892.0
  []
[]

[Functions]
  [q3_func]
    type = PiecewiseLinear
    x = '0.00000 0.02275 0.06825 0.11375 0.15925 0.20475 0.25025 0.29575 0.34125 0.38675 0.43225 0.45500'
    y = '2.0415e+07 2.0415e+07 2.8019e+07 3.4112e+07 3.8367e+07 4.0553e+07 4.0553e+07 3.8367e+07 3.4112e+07 2.8019e+07 2.0415e+07 2.0415e+07'
    axis = y
  []
  [T_bulk_func]
    type = PiecewiseLinear
    x = '0.00000 0.04550 0.09100 0.13650 0.18200 0.22750 0.27300 0.31850 0.36400 0.40950 0.45500'
    y = '892.000 898.711 914.631 935.054 958.879 984.821 1011.481 1037.423 1061.248 1081.671 1097.592'
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
    coefficient = 1161.28
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
  # Nb-1Zr Clad: k = 41.8 + 0.0125 * T [Touloukian 1970]
  [k_clad]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = 'T'
    expression = '41.8 + 0.0125 * T'
    block = '2'
  []
[]

[Postprocessors]
  [vol_fuel]
    type = VolumePostprocessor
    block = '1'
  []
  [vol_clad]
    type = VolumePostprocessor
    block = '2'
  []
  [Q_removed_boundary]
    type = SideDiffusiveFluxIntegral
    diffusivity = thermal_conductivity
    variable = T
    boundary = 'clad_right'
  []
  [T_fuel_max]
    type = NodalExtremeValue
    variable = T
    value_type = max
    block = '1'
  []
  [T_fuel_avg]
    type = ElementAverageValue
    variable = T
    block = '1'
  []
  [T_clad_max]
    type = NodalExtremeValue
    variable = T
    value_type = max
    block = '2'
  []
  [T_clad_avg]
    type = ElementAverageValue
    variable = T
    block = '2'
  []
  [T_fuel_exit_center]
    type = PointValue
    variable = T
    point = '0.0 0.45500 0.0'
  []
  [T_clad_exit_surf]
    type = PointValue
    variable = T
    point = '0.00633 0.45500 0.0'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
[]

[Outputs]
  csv = true
[]
