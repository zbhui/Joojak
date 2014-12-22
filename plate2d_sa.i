[GlobalParams]
 order = FIRST
 family = MONOMIAL
  	
  mach = 0.2
  reynolds = 2E+05
  	
  variables = 'rho momentum_x momentum_y momentum_z rhoe rhon'
[]

[Mesh]
  type = FileMesh
  file = ../high-order-workshop/C1.4_plate/a2-125-2s.msh
  dim = 2

  boundary_id = '1 2 3 4 5' 
  boundary_name = 'symmetric wall right top left'

  block_id = '0'
  block_name = 'fluid'
[]


[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    #petsc_options = '-ksp_monitor -ksp_view -snes_test_display'
    #petsc_options_iname = '-pc_type -snes_type'
    petsc_options_iname = '-ksp_type  -pc_type'
    petsc_options_value = 'gmres       ilu'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = newton
  num_steps = 1000
  l_tol = 1e-02
  l_max_its = 100
 	
  nl_max_its = 100
  nl_rel_tol = 1e-04

  [./TimeStepper]
    type = RatioTimeStepper
    dt = 100
    ratio = 2
    step = 2
    max_dt = 100	
  [../]
[]

[CFDVariables]
[]

[./CFDAuxVariables]
  type = SAAuxVariable
  aux_variables = 'pressure mach velocity_x velocity_y velocity_z eddy_viscosity'    
[../]

[CFDICs]
  type = SAIC
[]

[CFDKernels]
  type = SACellKernel
[]

[CFDDGKernels]
  type = SAFaceKernel
[]

[CFDBCs]
  type = SABC
  boundary = 'symmetric wall right top left'
[]

[CFDPostprocessor]
[]

[Outputs]
  [./exodus]
    type = Exodus
    output_initial = true
    interval = 1 					#间隔
  [../]
  [./console]
    type = Console	
    perf_log = true
    linear_residuals = true
    nonlinear_residuals =  true	
  [../]
[]

[Materials]
  [./cell_material]
  block = fluid
    type = SACellMaterial
    wall_distance = distance
  [../]

  [./face_material]
    block = fluid
    type = SAFaceMaterial
  [../]

  [./wall_material]
    boundary = wall
    bc_type = isothermal_wall
    type = SABndMaterial
  [../]
  [./far_field_material]
    boundary = 'left right top'
    bc_type = far_field
    type = SABndMaterial
  [../]
  [./symmetric_material]
    boundary = 'symmetric'
    bc_type = symmetric
    type = SABndMaterial
  [../]

[]

