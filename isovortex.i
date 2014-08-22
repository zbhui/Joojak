# 全局变量
[GlobalParams]
 	order = SECOND
 	family = MONOMIAL
  	
  gamma = 1.4
  mach = 0.1
  reynolds = 10.0
  prandtl = 0.72
  	
  attack = 0
  slide = 0
  	
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
	lumping = false
[]

# 网格
[Mesh]
  type = GeneratedMesh
  dim = 2
  
  nx = 20
  ny = 20
  
  xmin = -10
  xmax = 0

  ymin = -10
  ymax = 0
  
  block_id = '0'
  block_name = 'fluid'
[]

[AuxVariables]
  [./pressure]
  [../]

  [./mach]
  [../]

  [./velocity_x]
  [../]

  [./velocity_y]
  [../]

  [./velocity_z]
  [../]
[]

[AuxKernels]
  [./pressure]
		type = NSAuxVariable
		variable = pressure
  [../]

  [./mach]
		type = NSAuxVariable
		variable = mach
  [../]

  [./velocity_x]
		type = NSAuxVariable
		variable = velocity_x
  [../]

  [./velocity_y]
		type = NSAuxVariable
		variable = velocity_x
  [../]

  [./velocity_z]
		type = NSAuxVariable
		variable = velocity_x
  [../]
[]

[Preconditioning]
	[./SMP]
		type = SMP
		full = true

	  #petsc_options = '-ksp_monitor -ksp_view -snes_test_display'
    #petsc_options_iname = '-pc_type -snes_type'
  	#petsc_options_value = 'lu test'
    petsc_options_iname = '-pc_type '
  	petsc_options_value = 'ilu'
	[../]

[]
# 非线性系统求解
[Executioner]
  type = Transient
  solve_type = NEWTON
 	scheme = 'bdf2'
  dt = 0.001
  num_steps = 10
  
    # 线性迭代步的残差下降（相对）量级
 	l_tol = 1e-04
 # l_abs_step_tol = -1e-04
   # 最大线性迭代步	
 	l_max_its = 10
 	
 	# 最大非线性迭代步
 	nl_max_its = 10
 	# 非线性迭代的残值下降（相对）量级
  	nl_rel_tol = 1e-05
  	# 非线性迭代绝对残值
  	#nl_abs_tol = 1e-05

  	
	 abort_on_solve_fail = true	
  #end_time = 0.1
  
  	[./Adaptivity]
  	
 	[../]
[]

[Functions]
  [./exact_rho]
    type = IsoVortexExact
  [../]
[]

[Postprocessors]
  #[./h]
   # type = AverageElementSize
    #variable = rho
 # [../]

 # [./dofs]
  #  type = NumDOFs
 # [../]

  [./l2_err]
    type = ElementL2Error
    variable = rho
    function = exact_rho
  [../]
  
  #[./nodes]
  #  type = NumNodes
  #[../]

  #[./elements]
   # type = NumElems
  #[../]

 # [./residuals]
 #   type = Residual
 # [../]
  
 #[./integral_left]
 #  type = ElementIntegralVariablePostprocessor
 #  variable = rho
 #[../]  
[]

# 输出和后处理
[Outputs]
  	file_base = isovortex
	
 	
	[./exodus]
		type = Exodus
		output_initial = true
		
		interval = 1 					#间隔
		oversample = true
		refinements = 0
	[../]
	
	[./console]
		type = Console	
		perf_log = true
		linear_residuals = true
	  	nonlinear_residuals =  true	
		#verbose = true
    	#setup_log_early = true
    	#time_precision = 6
    	#fit_mode = 100
	[../]
	[./debug]
	    type = DebugOutput
  		#show_var_residual_norms = true
 		# show_actions = true
  		#show_top_residuals = 5
	[../]
[]

# 变量
[Variables]
	active = 'rho momentum_x momentum_y momentum_z rhoe'

	[./rho]
		[./InitialCondition] 
			type = IsoVortexIC
		[../]
  [../]

 	[./momentum_x]
		[./InitialCondition] 
			type = IsoVortexIC
		[../]
  [../]
  
 	[./momentum_y]
		[./InitialCondition] 
			type = IsoVortexIC
		[../]
  [../]
  	
   [./momentum_z]
		[./InitialCondition] 
			type = IsoVortexIC
		[../]
  [../]
  	
  [./rhoe]
		[./InitialCondition] 
			type = IsoVortexIC
		[../]
  [../]	
		
[]

# 体积分
[Kernels]
	[./mass_time]
		type = TimeDerivative
		variable = rho
	[../]		
	[./x-momentumum_time]
		type = TimeDerivative
		variable = momentum_x
	[../]	
	[./y-momentumum_time]
		type = TimeDerivative
		variable = momentum_y
	[../]
	[./z-momentumum_time]
		type = TimeDerivative
		variable = momentum_z
	[../]		
	[./total-energy_time]
		type = TimeDerivative
		variable = rhoe
	[../]

	[./mass_space]
		type = EulerCellKernel
		variable = rho
	[../]		
	[./x-momentumum_space]
		type = EulerCellKernel
		variable = momentum_x
	[../]	
	[./y-momentumum_space]
		type = EulerCellKernel
		variable = momentum_y
	[../]
	[./z-momentumum_space]
		type = EulerCellKernel
		variable = momentum_z
	[../]		
	[./total-energy_space]
		type = EulerCellKernel
		variable = rhoe
	[../]
[]


[DGKernels]
	[./mass_dg]
		type = EulerFaceKernel
		variable = rho
	[../]		
	[./x-momentumum_dg]
		type = EulerFaceKernel
		variable = momentum_x
	[../]	
	[./y-momentumum_dg]
		type = EulerFaceKernel
		variable = momentum_y
	[../]
	[./z-momentumum_dg]
		type = EulerFaceKernel
		variable = momentum_z
	[../]		
	[./total-energy_dg]
		type = EulerFaceKernel
		variable = rhoe
	[../]
[]

# 边界条件
[BCs]
	[./mass_bc]
		boundary = 'left right bottom top'
		type =EulerBC
		variable = rho
	[../]		
	[./x-momentumum_bc]
		boundary = 'left right bottom top'
		type =EulerBC
		variable = momentum_x
	[../]	
	[./y-momentumum_bc]
		boundary = 'left right bottom top'
		type =EulerBC
		variable = momentum_y
	[../]
	[./z-momentumum_bc]
		boundary = 'left right bottom top'
		type =EulerBC
		variable = momentum_z
	[../]		
	[./total-energy_bc]
		boundary = 'left right bottom top'
		type =EulerBC
		variable = rhoe
	[../]
[]

# 材料属性
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

