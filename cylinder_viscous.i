[Mesh]
  type = FileMesh
  file = grids/cylinder_fine.msh
  dim = 2
  
  block_id = 10
  block_name = 'fluid'
  
  boundary_id = '8 9'
  boundary_name = 'far_field wall'
[]

[MeshModifiers]
  [./extrude]
    type = BuildSideSetFromBlock
    new_boundary = boundary_from_block
  [../]
[]

[Problem]
  type = NavierStokesProblem
  mach = 0.1
  reynolds = 40
  aux_variables = distance
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

    [./partition]
      type = ProcessorIDAux 
      variables = proc_id
      order = CONSTANT
      family = MONOMIAL
    [../]

   [./distance]
      type = NearestNodeDistanceAux
      variables = distance
      order = FIRST
      family = LAGRANGE
      boundary = boundary_from_block
      paired_boundary = wall
    [../]

  [../]

[ICs]
  type = CFDPassFlowIC 
  velocity = 1
[]

[Adaptivity]
  marker = marker
  max_h_level = 2
  start_time = 200
  stop_time = 800
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
      coarsen = 0.1
      refine = 0.7
    [../]
  [../]
[]

[Materials]
  [./cell_material]
    block = 10
    type = CLawCellMaterial
    variables = 'rho momentum_x momentum_y momentum_z rhoe'
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
    bc_type = isothermal_wall
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

