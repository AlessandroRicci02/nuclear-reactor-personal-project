[Mesh]
  coord_type = RZ
  [cell]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 24
    ny = 160
    xmin = 0.010595
    xmax = 0.012895
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
    y = '1834460.26455 1922592.94367 2426880.56727 2687461.676 2672847.69042 2438370.2633 1926314.08136 1826312.18945 1826312.18945'
  []
  [inner_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '1023.15587875 1008.62444738 995.582383794 978.686184375 958.723277908 936.789828979 913.596931604 895.118123366 850'
  []
  [outer_bulk]
    type = PiecewiseLinear
    axis = y
    x = '0 0.0875 0.175 0.2625 0.35 0.4375 0.525 0.6125 0.7'
    y = '1023.90620417 1010.61016782 997.750018507 980.272051027 959.001543326 935.206499575 911.001502011 878.991448749 850'
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
    coefficient = 564.5151233670873
    T_infinity = inner_bulk
  []
  [outer]
    type = ConvectiveFluxFunction
    variable = T
    boundary = right
    coefficient = 645.0780629439561
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
