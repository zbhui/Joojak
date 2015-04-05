[GlobalParams]
  family = MONOMIAL
  mach = 0.2
  reynolds = 40
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
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
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
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
    output_initial = true
    interval = 1 					
  [../]
	
  [./console]
    type = Console	
    perf_log = true
    linear_residuals = true
    nonlinear_residuals =  true	
  [../]
[]

[AuxVariables]
  [./proc_id]
    order = FIRST
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./proc_id]
    type = ProcessorIDAux
    variable = proc_id
  [../]
[]

[CFDICs]
  type = CFDPassFlowIC
[]

[CFDBCs]
  type = NSBC
  boundary = '8 9' 
[]

[Postprocessors]
  [./force_total-x]
    type = CFDForcePostprocessor
    direction_by = x
    force_type = total
    boundary  = wall
  [../]
  [./force_total-y]
    type = CFDForcePostprocessor
    direction_by = y
    force_type = total
    boundary  = wall
  [../]
  [./force_total-z]
    type = CFDForcePostprocessor
    direction_by = z
    force_type = total
    boundary  = wall
  [../]
[]

[Materials]
  [./cell_material]
    block = 10
    type = NSCellMaterial
  [../]
  [./face_material]
    block = 10
    type = NSFaceMaterial
  [../]
  [./far_field_material]
    boundary = far_field
    bc_type = far_field
    type = NSBndMaterial
  [../]
  [./wall_material]
    boundary = wall
    bc_type = adiabatic_wall
    type = NSBndMaterial
  [../]
[]

