[Mesh]
  type = FileMesh
  file = grids/cylinder_tri.msh
  dim = 2
  
  block_id = 10
  block_name = 'fluid'
  
  boundary_id = '8 9'
  boundary_name = 'far_field wall'
[]

[MeshModifiers]
  [./fake_boudary]
    type = BuildSideSetFromBlock
  [../]
[]


[Problem]
  type = NavierStokesProblem
  order = FIRST
  family = MONOMIAL
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
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

  [./distance_to_left_nodes]
  [../]
  [./distance]
  [../]
[]


[AuxKernels]
  [./nodal_distance_aux]
    type = NearestNodeDistanceAux
    variable = distance_to_left_nodes
    boundary = boundary_from_block
    paired_boundary = 9
  [../]
  [./penetration_aux]
    type = PenetrationAux
    variable = distance
    boundary = boundary_from_block
    paired_boundary = 9
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
