[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
[]

[Problem]
  type = Riemann2DProblem
  aux_variables = artificial_vis
  [./Variables]
    order = SECOND
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
      marker = jump_marker
      order = FIRST
      family = MONOMIAL
    [../]
  [../]
[]

[ICs]
  type = CLawIC 
[]


[Adaptivity]
  initial_marker = marker
  marker = marker
  max_h_level = 2
  initial_steps = 2
  [./Indicators]
    [./error]
      type = FluxJumpIndicator
      variable = rhoe
    [../]
  [../]
  [./Markers]
    [./marker]
      type = ComboMarker
      markers = 'jump_marker box1 box2'
    [../]    
    [./jump_marker]
      type = ErrorFractionMarker
      indicator = error
      coarsen = 0.2
      refine = 0.8
    [../]
    [./box1]
      type = BoxMarker
      bottom_left = '0 0.4 0'
      top_right = '0.53 0.53 0'
      inside = refine
      outside = coarsen
    [../]
    [./box2]
      type = BoxMarker
      bottom_left = '0.47 0.0 0'
      top_right = '0.53 0.53 0'
      inside = refine
      outside = coarsen
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
  [./bc_material]
    type = CLawBoundaryMaterial
    boundary = ANY_BOUNDARY_ID
  [../]
[]

[Postprocessors]
[]

[UserObjects]
  [./layered_integral]
    type = CFDForceUserObject
    boundary = 0
    #execute_on = timestep_end
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
  num_steps = 1000
  #scheme = crank-nicolson
  scheme = bdf2
  l_tol = 1e-02
  #l_abs_step_tol = -1e-04
  l_max_its = 100
 	
  nl_max_its = 100
  nl_rel_tol = 1e-02
  #nl_abs_tol = 1e-03
  end_time = 0.2

  [./TimeStepper]
    type = RatioTimeStepper
    dt = 0.001
    ratio = 1
    step = 2
    max_dt = 100	
  [../]
[]

[Outputs]
  [./exodus]
    output_initial = true
    type = Exodus
    interval = 1 
    refinements = 0					
  [../]
	
  [./console]
    type = Console	
    perf_log = true
    output_on = 'linear nonlinear'
  [../]
[]

