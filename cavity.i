[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 40
  ny = 40
[]

[Problem]
  type = NavierStokesProblem
  mach = 0.1
  reynolds = 100

  [./Variables]
    order = FIRST
    family = MONOMIAL
    variables = 'rho momentum_x momentum_y momentum_z rhoe'
  [../]

  [./AuxVariables]
    [./Output]
      type = CFDAuxVariable 
      variables = 'pressure velocity_x velocity_y velocity_z mach'
      order = FIRST
      family = MONOMIAL
    [../]

  [../]
[]

[ICs]
  type = CFDPassFlowIC 
  velocity = 0
[]


[Materials]
  [./cell_material]
    block = ANY_BLOCK_ID
    type = CLawCellMaterial
    variables = 'rho momentum_x momentum_y momentum_z rhoe'
  [../]
  [./face_material]
    block = ANY_BLOCK_ID
    type = CLawFaceMaterial
  [../]
  [./far]
    type = CLawBoundaryMaterial
    boundary = top
    bc_type = far_field
  [../]
  [./wall]
    type = CLawBoundaryMaterial
    boundary = 'left right bottom'
    bc_type = isothermal_wall
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
    dt = 0.01
    ratio = 2
    step = 2
    max_dt = 100	
  [../]
[]

[Postprocessors]
[]

[Outputs]
  [./exodus]
    type = Exodus
    interval = 1 
    output_on = 'initial timestep_end'					
  [../]
	
  [./console]
    type = Console	
    perf_log = true
    output_on = 'linear nonlinear'
  [../]
[]



