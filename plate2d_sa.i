# 全局变量
[GlobalParams]
 	order = FIRST
 	family = MONOMIAL
  	
  mach = 0.2
  reynolds = 2E+05
  	
  variables = 'rho momentum_x momentum_y momentum_z rhoe rhon'
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
	[./potential]
		order = FIRST
		family = LAGRANGE
	[../]
	[./distance]
	[../]
[]

[AuxKernels]
  [./pressure]
		type = SAAuxVariable
		variable = pressure
  [../]

  [./mach]
		type = SAAuxVariable
		variable = mach
  [../]

  [./velocity_x]
		type = SAAuxVariable
		variable = velocity_x
  [../]

  [./velocity_y]
		type = SAAuxVariable
		variable = velocity_y
  [../]

  [./velocity_z]
		type = SAAuxVariable
		variable = velocity_z
  [../]
  [./eddy_viscosity]
		type = SAAuxVariable
		variable = eddy_viscosity
  [../]
	[./distace]
		type = NearestWallDistance
		variable = distance
		potential = potential
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
    petsc_options_iname = '-ksp_type  -pc_type -pc_sub_type'
  	petsc_options_value = 'gmres 				asm  ilu'
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
  	nl_rel_tol = 1e-02
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
[Dampers]
  [./u_damp]
    type = ConstantDamper
    variable = rhon
    damping = 1
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
	[../]
[]

[Debug]
[../]

# 变量
[Variables]
	[./rho]
		[./InitialCondition] 
			type = SAIC
		[../]
  [../]
 	[./momentum_x]
		[./InitialCondition] 
			type = SAIC
		[../]
  [../]  
 	[./momentum_y]
		[./InitialCondition] 
			type = SAIC
		[../]
  [../] 	
  [./momentum_z]
		[./InitialCondition] 
			type = SAIC
		[../]
  [../] 	
  [./rhoe]
		[./InitialCondition] 
			type = SAIC
		[../]
  [../]	
  [./rhon]
		[./InitialCondition] 
			type = SAIC
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
	[./rhon_time]
		type = TimeDerivative
		variable = rhon
	[../]

	[./mass_space]
		type = SACellKernel
		variable = rho
	[../]		
	[./x-momentumum_space]
		type = SACellKernel
		variable = momentum_x
	[../]	
	[./y-momentumum_space]
		type = SACellKernel
		variable = momentum_y
	[../]
	[./z-momentumum_space]
		type = SACellKernel
		variable = momentum_z
	[../]		
	[./total-energy_space]
		type = SACellKernel
		variable = rhoe
	[../]
	[./rhon_space]
		type = SACellKernel
		variable = rhon
	[../]
[]


[DGKernels]
	[./mass_dg]
		type = SAFaceKernel
		variable = rho
	[../]		
	[./x-momentumum_dg]
		type = SAFaceKernel
		variable = momentum_x
	[../]	
	[./y-momentumum_dg]
		type = SAFaceKernel
		variable = momentum_y
	[../]
	[./z-momentumum_dg]
		type = SAFaceKernel
		variable = momentum_z
	[../]		
	[./total-energy_dg]
		type = SAFaceKernel
		variable = rhoe
	[../]
	[./rhon_dg]
		type = SAFaceKernel
		variable = rhon
	[../]
[]

# 边界条件
[BCs]
	[./mass_bc]
		boundary = '1 2 3 4 5'
		type =SABC
		variable = rho
	[../]		
	[./x-momentumum_bc]
		boundary = '1 2 3 4 5'
		type =SABC
		variable = momentum_x
	[../]	
	[./y-momentumum_bc]
		boundary = '1 2 3 4 5'
		type =SABC
		variable = momentum_y
	[../]
	[./z-momentumum_bc]
		boundary = '1 2 3 4 5'
		type =SABC
		variable = momentum_z
	[../]		
	[./total-energy_bc]
		boundary = '1 2 3 4 5'
		type =SABC
		variable = rhoe
	[../]
	[./rhon_bc]
		boundary = '1 2 3 4 5'
		type =SABC
		variable = rhon
	[../]
[]

# 材料属性
[Materials]
  [./cell_material]
		block = fluid
    type = SACellMaterial
		wall_distance = distance
  [../]

  [./face_material]
		block = fluid
    type = SAFaceMaterial
  [../]

  [./wall_material]
		boundary = wall
		bc_type = isothermal_wall
    type = SABndMaterial
  [../]
  [./far_field_material]
		boundary = 'left right top'
		bc_type = far_field
    type = SABndMaterial
  [../]
  [./symmetric_material]
		boundary = 'symmetric'
		bc_type = symmetric
    type = SABndMaterial
  [../]

[]

