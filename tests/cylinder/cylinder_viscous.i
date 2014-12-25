[GlobalParams]
  order = FIRST
  family = MONOMIAL
  	
  mach = 0.1
  reynolds = 40	
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
[]

[Mesh]
  type = FileMesh
  file = grids/cylinder.msh
  dim = 2
  
  block_id = 10
  block_name = 'fluid'
  
  boundary_id = '8 9'
  boundary_name = 'far_field wall'
	
  uniform_refine = 0 
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
  num_steps = 1
  l_tol = 1e-02
  #l_abs_step_tol = -1e-04
  l_max_its = 100
 	
  nl_max_its = 100
  nl_rel_tol = 1e-04
  #nl_abs_tol = 1e-05

  [./TimeStepper]
    type = RatioTimeStepper
    dt = 0.100
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
    oversample = true
    refinements = 0
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

[CFDVariables]
[]

[./CFDAuxVariables]
  type = NSAuxVariable
  aux_variables = 'pressure mach velocity_x velocity_y velocity_z'    
[../]

[CFDICs]
  type = CFDPassFlowIC
[]

[CFDKernels]
  type = NSCellKernel
[]

[CFDDGKernels]
  type = NSFaceKernel
[]

[CFDBCs]
  type = NSBC
  boundary = '8 9' 
[]

[CFDPostprocessor]
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

