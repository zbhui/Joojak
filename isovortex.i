[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = -10
  xmax = 0
  ymin = -10
  ymax = 0
  block_id = '0'
  block_name = 'fluid'
[]

[Problem]
  type = EulerProblem
  order = FIRST
  family = MONOMIAL
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
  mach = 0.38
[]

[ICs]
  type = IsoVortexIC 
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
    variables = 'rho momentum_x momentum_y momentum_z rhoe'
  [../]
  [./face_material]
    block = 0
    type = CLawFaceMaterial
  [../]
  [./bc_material]
    type = IsoVortexBndMaterial
    boundary = '0 1 2 3'
  [../]
[]

[Postprocessors]
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

