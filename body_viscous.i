# 全局变量
[GlobalParams]
 	order = SECOND
 	family = MONOMIAL
  	
  mach = 0.5
  reynolds = 5000
	attack = 1
	sideslip = 0
	pitch = 0
	yaw = 180
	roll = 0
	ref_area = 0.05
  	
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
[]

# 网格
[Mesh]
  type = FileMesh
  file = ../high-order-workshop/C2.3_body/btc0-NLR-E1.v2.msh
  dim = 3

  boundary_id = '1 2 3'
  boundary_name = 'wall symmetric far_field'

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

    petsc_options_iname = 'ksp_type -pc_type '
  	petsc_options_value = 'bcgs bjacobi'
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
  	nl_rel_tol = 1e-3
  	# 非线性迭代绝对残值
  	nl_abs_tol = 1e-010

  	
	 #abort_on_solve_fail = true	
  #end_time = 0.1
  
	[./TimeStepper]
		type = RatioTimeStepper
		dt = 1E+02
		ratio = 2
		step = 2
		max_dt = 1E+02
	[../]
[]


[Postprocessors]
	[./res_evalutions]
		type = NumTimeStep
	[../]

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


	[./force_total-x]
  	type = CFDForcePostprocessor
		direction_by = x
		force_type = total
		boundary  = wall
	[../]
	[./force_total-y]
  	type = CFDForcePostprocessor
		direction_by = y
		force_type = total
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
	csv = true
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
		boundary = '1 2 3'
		type = NSBC
		variable = rho
	[../]		
	[./x-momentumum_bc]
		boundary = '1 2 3'
		type = NSBC
		variable = momentum_x
	[../]	
	[./y-momentumum_bc]
		boundary = '1 2 3 '
		type = NSBC
		variable = momentum_y
	[../]
	[./z-momentumum_bc]
		boundary = '1 2 3 '
		type = NSBC
		variable = momentum_z
	[../]		
	[./total-energy_bc]
		boundary = '1 2 3'
		type = NSBC
		variable = rhoe
	[../]
[]

# 材料属性
[Materials]
  [./cell_material]
		block = fluid
    type = NSCellMaterial
  [../]

  [./face_material]
		block = fluid
    type = NSFaceMaterial
  [../]

  [./wall_material]
		boundary = wall
		bc_type = isothermal_wall
    type = NSBndMaterial
  [../]
  [./far_field_material]
		boundary = '3'
		bc_type = far_field
    type = NSBndMaterial
  [../]
  [./symmetric_material]
		boundary = 2
		bc_type = symmetric
    type = NSBndMaterial
  [../]


[]

