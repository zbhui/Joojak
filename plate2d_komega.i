# 全局变量
[GlobalParams]
 	order = SECOND
 	family = MONOMIAL
  	
  mach = 0.2
  reynolds = 1E+05
  	
  variables = 'rho momentum_x momentum_y momentum_z rhoe rhok rhoo'
[]

# 网格
[Mesh]
  type = FileMesh
  file = ../high-order-workshop/C1.4_plate/a2-125-2s.msh
  dim = 2

  boundary_id = '1 2 3 4 5' 
  boundary_name = 'symmetric wall right top left'

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
  [./eddy_viscosity]
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
  [./eddy_viscosity]
		type = NSAuxVariable
		variable = eddy_viscosity
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
  	petsc_options_value = 'gmres 				lu  '
	[../]

[]
# 非线性系统求解
[Executioner]
  type = Transient
  solve_type = NEWTON
  num_steps = 100000
  
    # 线性迭代步的残差下降（相对）量级
 	l_tol = 1e-01
 # l_abs_step_tol = -1e-04
   # 最大线性迭代步	
 	l_max_its = 100
 	
 	# 最大非线性迭代步
 	nl_max_its = 10
 	# 非线性迭代的残值下降（相对）量级
  	nl_rel_tol = 1e-2
  	# 非线性迭代绝对残值
  	nl_abs_tol = 1e-010

  	#restart_file_base = plate2d_komega_checkpoint_cp/0026 
	 #abort_on_solve_fail = true	
  #end_time = 0.1
  
	[./TimeStepper]
		type = RatioTimeStepper
		dt = 1e-02
		ratio = 2
		step = 2
		max_dt = 1e+01
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
	[./force_form-z]
  	type = CFDForcePostprocessor
		direction_by = z
		force_type = form
		boundary  = wall
	[../]
	[./force_friction-z]
  	type = CFDForcePostprocessor
		direction_by = z
		force_type = friction
		boundary  = wall
	[../]
	[./force_total-z]
  	type = CFDForcePostprocessor
		direction_by = z
		force_type = total
		boundary  = wall
	[../]
 
[]

# 输出和后处理
[Outputs]
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
	[./rho]
		[./InitialCondition] 
			type = KOIC
		[../]
  [../]
 	[./momentum_x]
		[./InitialCondition] 
			type = KOIC
		[../]
  [../]  
 	[./momentum_y]
		[./InitialCondition] 
			type = KOIC
		[../]
  [../] 	
  [./momentum_z]
		[./InitialCondition] 
			type = KOIC
		[../]
  [../] 	
  [./rhoe]
		[./InitialCondition] 
			type = KOIC
		[../]
  [../]	
  [./rhok]
		[./InitialCondition] 
			type = KOIC
		[../]
  [../]	
  [./rhoo]
		[./InitialCondition] 
			type = KOIC
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
	[./rhok_time]
		type = TimeDerivative
		variable = rhok
	[../]
	[./rhoo_time]
		type = TimeDerivative
		variable = rhoo
	[../]

	[./mass_space]
		type = KOCellKernel
		variable = rho
	[../]		
	[./x-momentumum_space]
		type = KOCellKernel
		variable = momentum_x
	[../]	
	[./y-momentumum_space]
		type = KOCellKernel
		variable = momentum_y
	[../]
	[./z-momentumum_space]
		type = KOCellKernel
		variable = momentum_z
	[../]		
	[./total-energy_space]
		type = KOCellKernel
		variable = rhoe
	[../]
	[./rhok_space]
		type = KOCellKernel
		variable = rhok
	[../]
	[./rhoo_space]
		type = KOCellKernel
		variable = rhoo
	[../]
[]


[DGKernels]
	[./mass_dg]
		type = KOFaceKernel
		variable = rho
	[../]		
	[./x-momentumum_dg]
		type = KOFaceKernel
		variable = momentum_x
	[../]	
	[./y-momentumum_dg]
		type = KOFaceKernel
		variable = momentum_y
	[../]
	[./z-momentumum_dg]
		type = KOFaceKernel
		variable = momentum_z
	[../]		
	[./total-energy_dg]
		type = KOFaceKernel
		variable = rhoe
	[../]
	[./rhok_dg]
		type = KOFaceKernel
		variable = rhok
	[../]
	[./rhoo_dg]
		type = KOFaceKernel
		variable = rhoo
	[../]
[]

# 边界条件
[BCs]
	[./mass_bc]
		boundary = '1 2 3 4 5'
		type =KOBC
		variable = rho
	[../]		
	[./x-momentumum_bc]
		boundary = '1 2 3 4 5'
		type =KOBC
		variable = momentum_x
	[../]	
	[./y-momentumum_bc]
		boundary = '1 2 3 4 5'
		type =KOBC
		variable = momentum_y
	[../]
	[./z-momentumum_bc]
		boundary = '1 2 3 4 5'
		type =KOBC
		variable = momentum_z
	[../]		
	[./total-energy_bc]
		boundary = '1 2 3 4 5'
		type =KOBC
		variable = rhoe
	[../]
	[./rhok_bc]
		boundary = '1 2 3 4 5'
		type =KOBC
		variable = rhok
	[../]
	[./rhoo_bc]
		boundary = '1 2 3 4 5'
		type =KOBC
		variable = rhoo
	[../]
[]

# 材料属性
[Materials]
  [./cell_material]
		block = fluid
    type = KOCellMaterial
  [../]

  [./face_material]
		block = fluid
    type = KOFaceMaterial
  [../]

  [./wall_material]
		boundary = wall
		bc_type = isothermal_wall
		#bc_type = far_field
    type = KOBndMaterial
  [../]
  [./far_field_material]
		boundary = 'left right'
		bc_type = far_field
    type = KOBndMaterial
  [../]
  [./symmetric_material]
		boundary = 'symmetric top'
		bc_type = symmetric
    type = KOBndMaterial
  [../]

[]

