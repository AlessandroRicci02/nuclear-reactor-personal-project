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
    y = '1085225.07735 1150459.36093 1485139.51095 1654183.81232 1654745.61419 1482515.30532 1169800.37761 1085524.11589 1085524.11589'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0825 0.165 0.2475 0.33 0.4125 0.495 0.5775 0.66'
    y = '888.137769061 885.332862943 882.358362331 878.337820095 873.568440376 868.38337824 863.161194234 858.255844288 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0825 0.165 0.2475 0.33 0.4125 0.495 0.5775 0.66'
    y = '888.137769061 885.332862943 882.358362331 878.337820095 873.568440376 868.38337824 863.161194234 858.255844288 850'
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
