[Mesh]
  type = FileMesh
  file = cell.msh
[]
[Variables/T]
  initial_condition = 0
[]
[Kernels]
  [diffusion]
    type = HeatConduction
    variable = T
  []
  [mod_source]
    type = BodyForce
    variable = T
    value = 4135.387213900906
    block = '1'
  []
  [liner_source]
    type = BodyForce
    variable = T
    value = 0
    block = '2'
  []
[]
[Materials]
  [moderator]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 10.0
    block = 1
  []
  [liner]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 6.0
    block = 2
  []
[]
[BCs]
  [main]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 10
    coefficient = 1164.170665594267
    T_infinity = 0.0
  []
  [aux]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 11
    coefficient = 443.0982722752106
    T_infinity = 0
  []
[]
[Postprocessors]
  [Q_main]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    boundary = 10
    T_fluid = 0.0
    htc = 1164.170665594267
  []
  [Q_aux]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    boundary = 11
    T_fluid = 0
    htc = 443.0982722752106
  []
  [T_mod_avg]
    type = ElementAverageValue
    variable = T
    block = 1
  []
  [T_liner_avg]
    type = ElementAverageValue
    variable = T
    block = 2
  []
[]
[VectorPostprocessors/nodes]
  type = NodalValueSampler
  variable = T
  sort_by = id
[]
[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  l_tol = 1e-10
[]
[Outputs]
  csv = true
[]
