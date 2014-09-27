# 全局变量
[GlobalParams]
 	order = SECOND
 	family = MONOMIAL
  	
  mach = 0.8
	reynolds = 10
	attack = 1.25
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
[]

# 网格
[Mesh]
  type = FileMesh
  file = ../high-order-workshop/C1.3_naca0012/N0012-fine-quad.msh
  dim = 2
  
  block_id = 0
  block_name = 'fluid'
  
  boundary_id = '1 4 2 3'
  boundary_name = 'far_top far_bottom wall_top wall_bottom'
	
	uniform_refine = 0 
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
		variable = velocity_y
  [../]

  [./velocity_z]
		type = NSAuxVariable
		variable = velocity_z
  [../]
[]

[Preconditioning]
	[./SMP]
		type = SMP
		full = true

	  #petsc_options = '-ksp_monitor -ksp_view -snes_test_display'
    #petsc_options_iname = '-pc_type -snes_type'
  	#petsc_options_value = 'lu test'
		#petsc_options = '-pc_sor_symmetric'
    petsc_options_iname = '-ksp_type  -pc_type'
  	petsc_options_value = 'gmres 				bjacobi'
	[../]

[]
# 非线性系统求解
[Executioner]
  type = Transient
  solve_type = newton
 	#scheme = 'bdf2'
  num_steps = 1000
  
    # 线性迭代步的残差下降（相对）量级
 	l_tol = 1e-02
 # l_abs_step_tol = -1e-04
   # 最大线性迭代步	
 	l_max_its = 50
 	
 	# 最大非线性迭代步
 	nl_max_its = 10
 	# 非线性迭代的残值下降（相对）量级
  	nl_rel_tol = 1e-03
  	# 非线性迭代绝对残值
  	nl_abs_tol = 1e-012

  	
	 #abort_on_solve_fail = true	
  #end_time = 0.1
  
	[./TimeStepper]
		type = RatioTimeStepper
		dt = 0.0100000
		ratio = 2
		step = 2
		max_dt = 1000000
	[../]
[]

[Functions]
  [./exact_rho]
    type = IsoVortexExact
  [../]
[]

[Postprocessors]
	[./residual_final]
  	type = Residual
	[../]
  
	[./residual_initial]
  	type = CFDResidual
	[../]

	[./run_time]
  	type = RunTime
		time_type = alive
	[../]

	[./force_form-x]
  	type = CFDForcePostprocessor
		direction_by = x
		force_type = form
		boundary  = '2 3'
	[../]
	[./force_form-y]
  	type = CFDForcePostprocessor
		direction_by = y
		force_type = form
		boundary  = '2 3'
	[../]
	[./force_form-z]
  	type = CFDForcePostprocessor
		direction_by = z
		force_type = form
		boundary  = '2 3'
	[../]
[]

# 输出和后处理
[Outputs]
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
			type = CFDPassFlowIC
		[../]
  [../]

 	[./momentum_x]
		[./InitialCondition] 
			type = CFDPassFlowIC
		[../]
  [../]
  
 	[./momentum_y]
		[./InitialCondition] 
			type = CFDPassFlowIC
		[../]
  [../]
  	
   [./momentum_z]
		[./InitialCondition] 
			type = CFDPassFlowIC
		[../]
  [../]
  	
  [./rhoe]
		[./InitialCondition] 
			type = CFDPassFlowIC
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
		boundary = '1 2 3 4'
		type =EulerBC
		variable = rho
	[../]		
	[./x-momentumum_bc]
		boundary = '1 2 3 4'
		type =EulerBC
		variable = momentum_x
	[../]	
	[./y-momentumum_bc]
		boundary = '1 2 3 4'
		type =EulerBC
		variable = momentum_y
	[../]
	[./z-momentumum_bc]
		boundary = '1 2 3 4'
		type =EulerBC
		variable = momentum_z
	[../]		
	[./total-energy_bc]
		boundary = '1 2 3 4'
		type =EulerBC
		variable = rhoe
	[../]
[]

# 材料属性
[Materials]
  [./cell_material]
		block = fluid
    type = EulerCellMaterial
  [../]

  [./face_material]
		block = fluid
    type = EulerFaceMaterial
  [../]

 	[./far_field_material]
		boundary = '1 4'
		bc_type = far_field
    type = EulerBndMaterial
  [../]

  [./wall_material]
		boundary = '2 3'
		bc_type = wall
    type = EulerBndMaterial
  [../]

[]

