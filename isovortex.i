[GlobalParams]
  order = SECOND
  family = MONOMIAL
  	
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
  use_displaced_mesh = true
[]

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
  #second_order = true
  uniform_refine = 0
  displacements = 'disp_x disp_y'
[]


[AuxVariables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxKernels]
  [./disp_x]
    type = FunctionAux
    function = func_disp_x
    variable = disp_x
  [../]
  [./disp_y]
    type = FunctionAux
    function = func_disp_y
    variable = disp_y
  [../]
[]
[Functions]
  [./func_disp_x]
    type = ParsedFunction
    value = t
  [../]
  [./func_disp_y]
    type = ParsedFunction
    value = t
  [../]
[]

[Problem]
  kernel_coverage_check = false
[]

[CFDVariables]
[]

[./CFDAuxVariables]
  type = NSAuxVariable
  aux_variables = 'pressure mach velocity_x velocity_y velocity_z'    
[../]

[CFDICs]
  type = IsoVortexIC
[]

[CFDKernels]
  type = EulerCellKernel
[]

[CFDDGKernels]
  type = EulerFaceKernel
[]

[CFDBCs]
  type = EulerBC
    boundary = 'left right bottom top'
[]

[CFDPostprocessor]
  time_type = alive
[]


[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    #petsc_options = '-help'
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type'
    petsc_options_value = 'gmres  hypre parasails'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'crank-nicolson'
  dt = 0.02
  end_time = 1
  num_steps = 5
  l_tol = 1e-04
  l_max_its = 10
  nl_max_its = 10
  nl_rel_tol = 1e-08
[]

[Functions]
  [./exact_rho]
    type = IsoVortexExact
  [../]
[]


[Postprocessors]
  alive_time = true
  [./l2_err]
    type = ElementL2Error
    variable = rho
    function = exact_rho
  [../] 
[]

[Outputs]
  [./exodus]
    type = Exodus
    output_on = 'timestep_begin'
    interval = 1 				
    refinements = 0
  [../]
  [./console]
    type = Console	
    perf_log = true
    output_on = linear
  [../]
[]

[Materials]
  [./cell_materical]
    block = 0
    type = EulerCellMaterial
  [../]

  [./face_materical]
    block = 0
    type = EulerFaceMaterial
  [../]

  [./bnd_materical]
    boundary = 'left right bottom top'
    type = IsoVortexBndMaterial
  [../]
[]
