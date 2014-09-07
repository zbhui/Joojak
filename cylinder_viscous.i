# 全局变量
[GlobalParams]
 	order = THIRD
 	family = MONOMIAL
  	
  gamma = 1.4
  mach = 0.1
  reynolds = 40.0
  prandtl = 0.72
  	
  attack = 0
  slide = 0
  	
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
[]

# 网格
[Mesh]
  type = FileMesh
  file = grids/cylinder_fine.msh
  dim = 2
  
  block_id = 10
  block_name = 'fluid'
  
  boundary_id = '8 9'
  boundary_name = 'far_field wall'
	
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
    petsc_options_iname = '-ksp_type  -pc_type '
  	petsc_options_value = 'gmres 				bjacobi  '
	[../]

[]
# 非线性系统求解
[Executioner]
  type = Transient
  solve_type = NEWTON
  num_steps = 100000
  
    # 线性迭代步的残差下降（相对）量级
 	l_tol = 1e-02
 # l_abs_step_tol = -1e-04
   # 最大线性迭代步	
 	l_max_its = 100
 	
 	# 最大非线性迭代步
 	nl_max_its = 100
 	# 非线性迭代的残值下降（相对）量级
  	nl_rel_tol = 1e-4
  	# 非线性迭代绝对残值
  	nl_abs_tol = 1e-010

  	
	 #abort_on_solve_fail = true	
  #end_time = 0.1
  
	[./TimeStepper]
		type = RatioTimeStepper
		dt = 1E+08
		ratio = 2
		step = 2
		max_dt = 1e+08
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
		boundary  = wall
	[../]
	[./force_friction-x]
  	type = CFDForcePostprocessor
		direction_by = x
		force_type = friction
		boundary  = wall
	[../]
	[./force_total-x]
  	type = CFDForcePostprocessor
		direction_by = x
		force_type = total
		boundary  = wall
	[../]
	[./force_form-y]
  	type = CFDForcePostprocessor
		direction_by = y
		force_type = form
		boundary  = wall
	[../]
	[./force_friction-y]
  	type = CFDForcePostprocessor
		direction_by = y
		force_type = friction
		boundary  = wall
	[../]
	[./force_total-y]
  	type = CFDForcePostprocessor
		direction_by = y
		force_type = total
		boundary  = wall
	[../]
 
[]

# 输出和后处理
[Outputs]
	file_base = cylinder_viscous
	[./exodus]
		type = Exodus
		output_initial = true
		
		interval = 1 					#间隔
		oversample = true
		refinements = 1
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
	[./checkpoint]
		type  = Checkpoint
		interval = 1 					#间隔
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
		type = NSCellKernel
		variable = rho
	[../]		
	[./x-momentumum_space]
		type = NSCellKernel
		variable = momentum_x
	[../]	
	[./y-momentumum_space]
		type = NSCellKernel
		variable = momentum_y
	[../]
	[./z-momentumum_space]
		type = NSCellKernel
		variable = momentum_z
	[../]		
	[./total-energy_space]
		type = NSCellKernel
		variable = rhoe
	[../]
[]


[DGKernels]
	[./mass_dg]
		type = NSFaceKernel
		variable = rho
	[../]		
	[./x-momentumum_dg]
		type = NSFaceKernel
		variable = momentum_x
	[../]	
	[./y-momentumum_dg]
		type = NSFaceKernel
		variable = momentum_y
	[../]
	[./z-momentumum_dg]
		type = NSFaceKernel
		variable = momentum_z
	[../]		
	[./total-energy_dg]
		type = NSFaceKernel
		variable = rhoe
	[../]
[]

# 边界条件
[BCs]
	[./mass_bc]
		boundary = '8 9'
		type =NSBC
		variable = rho
	[../]		
	[./x-momentumum_bc]
		boundary = '8 9'
		type =NSBC
		variable = momentum_x
	[../]	
	[./y-momentumum_bc]
		boundary = '8 9'
		type =NSBC
		variable = momentum_y
	[../]
	[./z-momentumum_bc]
		boundary = '8 9'
		type =NSBC
		variable = momentum_z
	[../]		
	[./total-energy_bc]
		boundary = '8 9'
		type =NSBC
		variable = rhoe
	[../]
[]

# 材料属性
[Materials]
  [./cell_material]
		block = 10
    type = NSCellMaterial
  [../]

  [./face_material]
		block = 10
    type = NSFaceMaterial
  [../]

 	[./far_field_material]
		boundary = far_field
		bc_type = far_field
    type = NSBndMaterial
  [../]

  [./wall_material]
		boundary = wall
		bc_type = wall
    type = NSBndMaterial
  [../]

[]

