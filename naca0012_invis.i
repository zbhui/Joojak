[Mesh]
  type = FileMesh
  file = ./N0012-coarse-quad.msh
  dim = 2
  
  block_id = 0
  block_name = 'fluid'
  
  boundary_id = '1 4 2 3'
  boundary_name = 'far_top far_bottom wall_top wall_bottom'
[]

[Problem]
  type = EulerProblem
  mach = 0.85
  attack = 1.0
  aux_variables = artificial_vis
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
    [./artificial_vis]
      type = ArtificialViscosityAuxKernel 
      variables = artificial_vis
      indicator = error
      marker = marker
      order = FIRST
      family = MONOMIAL
    [../]
  [../]

[ICs]
  type = CFDPassFlowIC 
  velocity = 1
[]



[Adaptivity]
  [./Indicators]
    [./error]
      type = FluxJumpIndicator
      variable = rho
    [../]
  [../]
  [./Markers]
    [./marker]
      type = ErrorFractionMarker
      indicator = error
      coarsen = 0.7
      refine = 0.9
    [../]
  [../]
[]

[Materials]
  [./cell_material]
    block = 0
    type = CLawCellMaterial
    variables = 'rho momentum_x momentum_y momentum_z rhoe'
  [../]
  [./face_material]
    block = 0
    type = CLawFaceMaterial
  [../]
  [./far_field_material]
    boundary = 'far_top far_bottom'
    bc_type = far_field
    type = CLawBoundaryMaterial
  [../]
  [./wall_material]
    boundary ='wall_top wall_bottom'
    bc_type = wall
    type = CLawBoundaryMaterial
  [../]
[]

[Postprocessors]
[]


[UserObjects]
  [./cfd_force]
    type = CFDForceUserObject
    boundary = 'wall_top wall_bottom'
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
    dt = 0.0100
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

