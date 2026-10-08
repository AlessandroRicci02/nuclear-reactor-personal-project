# MOOSE / Cardinal 2D Axisymmetric FEM YH2.0 Moderator Model: Peak Moderator Cell
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
    ymax = 0.46000
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
    initial_condition = 880.0
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
    value = 8.8140e+05
    block = '1'
  []
  # Distributed auxiliary cooling sink modeling 6 aux channels (880 K)
  [aux_cooling_sink]
    type = CoefReaction
    variable = T
    coefficient = 317418.18
    block = '1'
  []
  [aux_cooling_ref]
    type = BodyForce
    variable = T
    value = 279327995.32
    block = '1'
  []
[]

[BCs]
  # Convection and liner conduction from hot fuel annulus (T_bulk ~ 1000 K avg)
  [inner_convection]
    type = ConvectiveFluxFunction
    variable = T
    boundary = 'mod_left'
    T_infinity = 1000.0
    coefficient = 1373.11
  []
  # Symmetry boundary at unit cell boundary (mod_right)
[]

[Materials]
  # YH2.0 Moderator: k ~ 17.5 W/(m K) [Flint 1965, Mueller 1968, Hu 2021]
  [k_yh2]
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
