[GlobalParams]
  use_displaced_mesh = true
[]

[Mesh]
  file = ./grids/cylinder_fine.msh
  dim = 2
  
  block_id = 10
  block_name = 'fluid'
  
  boundary_id = '8 9'
  boundary_name = 'far_field wall'
  displacements = 'disp_x disp_y'
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]

  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  [./solid_x]
    type = SolidMechX
    variable = disp_x
    y = disp_y
  [../]

  [./solid_y]
    type = SolidMechY
    variable = disp_y
    x = disp_x
  [../]
[]

[Functions]
  [./func_x]
    type = ParsedFunction
    value = 1.6*y
  [../]
  [./func_y]
    type = ParsedFunction
    value = 0
  [../]
[]

[BCs]
  [./wall_x]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = wall
    function = func_x
  [../]

  [./wall_y]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = wall
    function = func_y
  [../]


[]

[Materials]
  [./constant]
    type = LinearElasticityMaterial
    block = fluid
    youngs_modulus = 1
    poissons_ratio = 0.3
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type  -pc_type'
    petsc_options_value = 'gmres       lu'
  [../]
[]

[Executioner]
  type = Steady

  #Preconditioned JFNK (default)
  solve_type = 'newton'
  l_tol = 1e-03
  l_max_its = 100
  nl_max_its = 100
  nl_rel_tol = 1e-08
[]

[Outputs]
  file_base = out
  interval = 1
  exodus = true
  output_on = 'initial timestep_end'
  [./console]
    type = Console
    perf_log = true
    output_on = 'timestep_end failed nonlinear linear'
  [../]
[]
