[Mesh]
  type = FileMesh
  file = ../high-order-workshop/C1.4_plate/a2-125-2s.msh
  dim = 2

  boundary_id = '1 2 3 4 5' 
  boundary_name = 'symmetric wall right top left'

  block_id = '0'
  block_name = 'fluid'
[]

[Mesh]
  type = FileMesh
  file = grids/cylinder_fine.msh
  dim = 2
  
  block_id = 10
  block_name = 'fluid'
  
  boundary_id = '8 9'
  boundary_name = 'far_field wall'
[]

[Problem]
  type = NavierStokesProblem
  order = FIRST
  family = MONOMIAL
  variables = 'rho momentum_x momentum_y momentum_z rhoe rhon'
  mach = 0.2
  reynolds = 40
[]

[ICs]
  type = CFDPassFlowIC 
  velocity = 1
[]

[AuxVariables]
  type = NSAuxVariable 
  aux_variables = 'pressure velocity_x velocity_y velocity_z mach'
  order = FIRST
  family = MONOMIAL
[]

[Materials]
  [./cell_material]
    block = 10
    type = CLawCellMaterial
    variables = 'rho momentum_x momentum_y momentum_z rhoe rhon'
  [../]
  [./face_material]
    block = 10
    type = CLawFaceMaterial
  [../]
  [./far_field_material]
    boundary = far_field
    bc_type = far_field
    type = CLawBoundaryMaterial
  [../]
  [./wall_material]
    boundary = wall
    bc_type = adiabatic_wall
    type = CLawBoundaryMaterial
  [../]
[]

[Postprocessors]
[]


[UserObjects]
  [./cfd_force]
    type = CFDForceUserObject
    boundary = wall
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    #petsc_options = '-ksp_monitor -ksp_view -snes_test_display'
    #petsc_options_iname = '-pc_type -snes_type'
    petsc_options_iname = '-ksp_type  -pc_type'
    petsc_options_value = 'gmres       lu'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = newton
  num_steps = 1000
  l_tol = 1e-02
  #l_abs_step_tol = -1e-04
  l_max_its = 100
 	
  nl_max_its = 100
  nl_rel_tol = 1e-04
  #nl_abs_tol = 1e-05

  [./TimeStepper]
    type = RatioTimeStepper
    dt = 100
    ratio = 2
    step = 2
    max_dt = 100	
  [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    interval = 1 					
  [../]
	
  [./console]
    type = Console	
    perf_log = true
    output_on = 'linear nonlinear'
  [../]
[]

[Problem]
  type = SAProblem
  order = FIRST
  family = MONOMIAL
  variables = 'rho momentum_x momentum_y momentum_z rhoe rhon'
  mach = 0.2
  reynolds = 1E+06
[]

[ICs]
  type = CFDPassFlowIC 
  velocity = 1
[]

[AuxVariables]
  type = NSAuxVariable 
  aux_variables = 'pressure velocity_x velocity_y velocity_z mach'
  order = FIRST
  family = MONOMIAL
[]

[Materials]
  [./cell_material]
    block = 0
    type = CLawCellMaterial
    variables = 'rho momentum_x momentum_y momentum_z rhoe rhon'
  [../]
  [./face_material]
    block = 0
    type = CLawFaceMaterial
  [../]
  [./far_field_material]
    boundary = far_field
    bc_type = far_field
    type = CLawBoundaryMaterial
  [../]
  [./wall_material]
    boundary = wall
    bc_type = adiabatic_wall
    type = CLawBoundaryMaterial
  [../]
[]

[Postprocessors]
[]


[UserObjects]
  [./cfd_force]
    type = CFDForceUserObject
    boundary = wall
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    #petsc_options = '-ksp_monitor -ksp_view -snes_test_display'
    #petsc_options_iname = '-pc_type -snes_type'
    petsc_options_iname = '-ksp_type  -pc_type'
    petsc_options_value = 'gmres       lu'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = newton
  num_steps = 1000
  l_tol = 1e-02
  #l_abs_step_tol = -1e-04
  l_max_its = 100
 	
  nl_max_its = 100
  nl_rel_tol = 1e-04
  #nl_abs_tol = 1e-05

  [./TimeStepper]
    type = RatioTimeStepper
    dt = 100
    ratio = 2
    step = 2
    max_dt = 100	
  [../]
[]

[Outputs]
  [./exodus]
    type = Exodus
    interval = 1 					
  [../]
	
  [./console]
    type = Console	
    perf_log = true
    output_on = 'linear nonlinear'
  [../]
[]



