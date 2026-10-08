[Mesh]
  type = FileMesh
  file = /home/aleri/OpenMC_projects/first_test/new_possible_best_configurations/automated_search/followup/thermal/cell_yh2.msh
[]
[Variables]
  [T]
    initial_condition = 1000
  []
[]
[Functions]
  [heat]
    type = ConstantFunction
    value = 1000000.0
  []
  [main_bulk]
    type = ConstantFunction
    value = 1000
  []
  [aux_bulk]
    type = ConstantFunction
    value = 1000
  []
[]
[Kernels]
  [conduction]
    type = HeatConduction
    variable = T
  []
  [source]
    type = BodyForce
    variable = T
    function = heat
  []
[]
[BCs]
  [main]
    type = ConvectiveFluxFunction
    variable = T
    boundary = '2'
    coefficient = 576.0895263321759
    T_infinity = main_bulk
  []
  [aux]
    type = ConvectiveFluxFunction
    variable = T
    boundary = '3'
    coefficient = 332.08913425672966
    T_infinity = aux_bulk
  []
[]
[Materials]
  [conductivity]
    type = GenericConstantMaterial
    prop_names = thermal_conductivity
    prop_values = 10.0
  []
  [main_bulk]
    type = GenericConstantMaterial
    prop_names = main_bulk
    prop_values = 1000
  []
  [aux_bulk]
    type = GenericConstantMaterial
    prop_names = aux_bulk
    prop_values = 1000
  []
[]
[Postprocessors]
  [T_YH2_max_K]
    type = NodalExtremeValue
    variable = T
    value_type = max
  []
  [T_YH2_min_K]
    type = NodalExtremeValue
    variable = T
    value_type = min
  []
  [Q_main_W_per_m]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    T_fluid = main_bulk
    htc = 576.0895263321759
    boundary = '2'
  []
  [Q_aux_W_per_m]
    type = ConvectiveHeatTransferSideIntegral
    T_solid = T
    T_fluid = aux_bulk
    htc = 332.08913425672966
    boundary = '3'
  []
[]
[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-9
[]
[Outputs]
  csv = true
  exodus = true
[]
