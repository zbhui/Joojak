[GlobalParams]
  order = FIRST
  family = MONOMIAL
  	
  mach = 0.1
  reynolds = 10	
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
[]

# 网格
[Mesh]
  type = GeneratedMesh
  dim = 2
  
  nx = 10
  ny = 5
  
  xmin = 0
  xmax = 4

  ymin = 0
  ymax = 2
  
  block_id = '0'
  block_name = 'fluid'
[]

[Functions]
  [./exact_rho]
    type = CouetteFlowExact
  [../]
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
  num_steps = 1
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

[Postprocessors]
  [./l2_err]
    type = ElementL2Error
    variable = rho
    function = exact_rho
  [../]

  [./residuals]
    type = Residual
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
  type = NSCellKernel
[]

[CFDDGKernels]
  type = NSFaceKernel
[]

[CFDBCs]
  type = NSBC
  boundary = '0 1 2 3' 
[]

[CFDPostprocessor]
[]

[Materials]
  [./cell_materical]
    block = 0
    type = NSCellMaterial
  [../]
  [./face_materical]
    block = 0
    type = NSFaceMaterial
  [../]
  [./bnd_materical]
    boundary = 'left right bottom top'
    type = CouetteFlowBndMaterial
  [../]
[]

