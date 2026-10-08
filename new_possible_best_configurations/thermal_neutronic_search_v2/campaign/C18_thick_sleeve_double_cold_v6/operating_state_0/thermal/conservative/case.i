[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010195000000000001
    xmax = 0.012995000000000001
    ymin = 0
    ymax = 0.7
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
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '1600317.89135 1566370.06808 1961744.38787 2183089.45846 2178310.44824 1971149.67283 1565818.71811 1606696.5523 1606696.5523'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '983.806322472 975.42499786 968.235173738 959.325302563 949.185105857 938.521220482 928.125076393 919.252090945 900'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '985.641141085 978.063973537 971.230221994 962.052112886 950.982914261 938.811342685 926.127550942 912.723411775 900'
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
    coefficient = 564.0290351222121
    T_infinity = inner_bulk
  []
  [outer]
    type = ConvectiveFluxFunction
    variable = T
    boundary = right
    coefficient = 644.9693792092712
    T_infinity = outer_bulk
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
