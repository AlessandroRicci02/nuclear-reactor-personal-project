[Mesh]
  coord_type = RZ
  [fuel]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 60
    xmin = 0
    xmax = 0.00585
    ymin = 0
    ymax = 0.455
    boundary_name_prefix = fuel
  []
  [fuel_id]
    type = SubdomainIDGenerator
    input = fuel
    subdomain_id = 1
  []
  [clad]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 8
    ny = 60
    xmin = 0.00585
    xmax = 0.0063300000000000006
    ymin = 0
    ymax = 0.455
    boundary_name_prefix = clad
  []
  [clad_id]
    type = SubdomainIDGenerator
    input = clad
    subdomain_id = 2
  []
  [stitch]
    type = StitchMeshGenerator
    inputs = 'fuel_id clad_id'
    stitch_boundaries_pairs = 'fuel_right clad_left'
  []
[]
[Variables/T]
  initial_condition = 1000
[]
[Functions]
  [qfuel]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '0 0.0455 0.091 0.1365 0.182 0.2275 0.273 0.3185 0.364 0.4095 0.455'
    y = '34751689.229096 34751689.229096 34751689.229096 34751689.229096 34751689.229096 34751689.229096 34751689.229096 34751689.229096 34751689.229096 34751689.229096 34751689.229096'
  []
  [qclad]
    type = PiecewiseConstant
    axis = y
    direction = left
    x = '0 0.0455 0.091 0.1365 0.182 0.2275 0.273 0.3185 0.364 0.4095 0.455'
    y = '1196603.32416 1196603.32416 1196603.32416 1196603.32416 1196603.32416 1196603.32416 1196603.32416 1196603.32416 1196603.32416 1196603.32416 1196603.32416'
  []
  [bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0455 0.091 0.1365 0.182 0.2275 0.273 0.3185 0.364 0.4095 0.455'
    y = '892 915.86005351059 938.25315740165 959.83976266097 980.98297161387 1001.8824104266 1022.6478294461 1043.3395671488 1063.9907963838 1084.6197548991 1105.2364694323'
  []
[]
[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
  []
  [fuel_source]
    type = BodyForce
    variable = T
    block = 1
    function = qfuel
  []
  [clad_source]
    type = BodyForce
    variable = T
    block = 2
    function = qclad
  []
[]
[BCs/convection]
  type = ConvectiveFluxFunction
  variable = T
  boundary = clad_right
  coefficient = 1164.1706655942671
  T_infinity = bulk
[]
[Materials]
  [bulk_temperature]
    type = GenericFunctionMaterial
    prop_names = bulk_temperature
    prop_values = bulk
    block = 2
  []
  [fuel]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = T
    expression = '1.37*T^0.41'
    block = 1
  []
  [clad]
    type = ParsedMaterial
    property_name = thermal_conductivity
    coupled_variables = T
    expression = '41.8+0.0125*T'
    block = 2
  []
[]
[Postprocessors]
  [T_fuel_max]
    type = NodalExtremeValue
    variable = T
    value_type = max
    block = 1
  []
  [T_clad_max]
    type = NodalExtremeValue
    variable = T
    value_type = max
    block = 2
  []
  [Q_removed]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    T_fluid = bulk_temperature
    htc = 1164.1706655942671
    boundary = clad_right
  []
[]
[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8
[]
[Outputs]
  csv = true
[]
