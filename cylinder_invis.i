[GlobalParams]
  order = FIRST
  family = MONOMIAL
  	
  mach = 0.38
  attack = 90
  	
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
[]

[Mesh]
  type = FileMesh
  file = grids/cylinder_quad.msh
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
    petsc_options_value = 'gmres       bjacobi'
  [../]

[]

[Executioner]
  type = Transient
  solve_type = newton
  #scheme = 'bdf2'
  num_steps = 1000
  l_tol = 1e-01
  #l_abs_step_tol = -1e-04
  l_max_its = 50
 	
  nl_max_its = 5
  nl_rel_tol = 1e-02
  #nl_abs_tol = 1e-05

  [./TimeStepper]
    type = RatioTimeStepper
    dt = 0.01
    ratio = 2
    step = 2
    max_dt = 1	
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
  type = EulerCellKernel
[]

[CFDDGKernels]
  type = EulerFaceKernel
[]

[CFDBCs]
  type = EulerBC
  boundary = '8 9' 
[]

[CFDPostprocessor]
[]

[Materials]
  [./cell_material]
    block = 10
    type = EulerCellMaterial
  [../]
  [./face_material]
    block = 10
    type = EulerFaceMaterial
  [../]
  [./far_field_material]
    boundary = far_field
    bc_type = far_field
    type = EulerBndMaterial
  [../]
  [./wall_material]
    boundary = wall
    bc_type = wall
    type = EulerBndMaterial
  [../]
[]

