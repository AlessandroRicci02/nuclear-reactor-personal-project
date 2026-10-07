# MOOSE 2D Axisymmetric FEM ZrH1.2 Moderator Model: AVG_CELL_MODERATOR
[Mesh]
  coord_type = RZ
  [mod_mesh]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 20
    ny = 60
    xmin = 0.00727
    xmax = 0.01197
    ymin = 0.0
    ymax = 0.45500
    boundary_name_prefix = mod
  []
  [mod_id]
    type = SubdomainIDGenerator
    input = mod_mesh
    subdomain_id = 1
  []
[]

[Variables]
  [T]
    initial_condition = 800.0
  []
[]

[Functions]
  [q3_mod_func]
    type = PiecewiseLinear
    x = '0.00000 0.02275 0.06825 0.11375 0.15925 0.20475 0.25025 0.29575 0.34125 0.38675 0.43225 0.45500'
    y = '1.8026e+04 1.8026e+04 2.4740e+04 3.0121e+04 3.3878e+04 3.5808e+04 3.5808e+04 3.3878e+04 3.0121e+04 2.4740e+04 1.8026e+04 1.8026e+04'
    axis = y
  []
  [T_bulk_func]
    type = PiecewiseLinear
    x = '0.00000 0.04550 0.09100 0.13650 0.18200 0.22750 0.27300 0.31850 0.36400 0.40950 0.45500'
    y = '800.000 811.063 837.309 870.979 910.256 953.024 996.976 1039.744 1079.021 1112.691 1138.937'
    axis = y
  []
[]

[Kernels]
  [heat_conduction]
    type = HeatConduction
    variable = T
  []
  [nuclear_heating]
    type = BodyForce
    variable = T
    function = q3_mod_func
    block = '1'
  []
  # Distributed auxiliary cooling sink modeling 6 aux channels
  [aux_cooling_sink]
    type = CoefReaction
    variable = T
    coefficient = 169157.02
    block = '1'
  []
  [aux_cooling_ref]
    type = BodyForce
    variable = T
    function = T_bulk_func
    value = 169157.02
    block = '1'
  []
[]

[BCs]
  # Convection at the inner liner boundary to the annular coolant channel
  [inner_convection]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'mod_left'
    T_infinity = T_bulk_func
    coefficient = 1141.19
  []
  # Symmetry boundary at cell boundary (mod_right)
[]

[Materials]
  # ZrH1.2 Moderator: k ~ 17.5 W/(m K)
  [k_zrh]
    type = GenericConstantMaterial
    prop_names = 'thermal_conductivity'
    prop_values = '17.50'
    block = '1'
  []
[]

[Postprocessors]
  [vol_mod]
    type = VolumePostprocessor
    block = '1'
  []
  [T_mod_max]
    type = NodalExtremeValue
    variable = T
    value_type = max
    block = '1'
  []
  [T_mod_avg]
    type = ElementAverageValue
    variable = T
    block = '1'
  []
  [T_mod_min]
    type = NodalExtremeValue
    variable = T
    value_type = min
    block = '1'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_type'
  petsc_options_value = 'hypre boomeramg gmres'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-8
[]

[Outputs]
  csv = true
[]
