# 全局变量
[GlobalParams]

  	
[]

[CFD]
  order = FIRST
  family = MONOMIAL
  type = EULER
  init_cond = IsoVortexIC
  #aux_variables = 'pressure mach velocity_x velocity_y velocity_z'
  boundary = 'left right bottom top'
  [./Material]
    [./cell_materical1]
      block = 0
      type = EulerCellMaterial
      variables = 'rho momentum_x momentum_y momentum_z rhoe'
    [../]
    [./face_materical]
      block = 0
      type = EulerFaceMaterial
      variables = 'rho momentum_x momentum_y momentum_z rhoe'
    [../]
    [./bnd_materical]
      boundary = 'left right bottom top'
      type = IsoVortexBndMaterial
      variables = 'rho momentum_x momentum_y momentum_z rhoe'
    [../]
  [../]
[]

# 网格
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
  uniform_refine = 0
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
  dt = 0.002
  end_time = 1
  num_steps = 1
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
  [./l2_err]
    type = ElementL2Error
    variable = rho
    function = exact_rho
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









