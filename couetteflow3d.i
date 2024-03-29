# 全局变量
[GlobalParams]
 	order = THIRD
 	family = MONOMIAL
  	
  mach = 0.1
  reynolds = 10.0
  	
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
[]

# 网格
[Mesh]
  type = GeneratedMesh
  dim = 3
  
  nx = 10
  ny = 5
	nz = 5
  
  xmin = 0
  xmax = 4

  ymin = 0
  ymax = 2

  zmin = 0
  zmax = 2
  
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
  dt = 1E+10
  num_steps = 1000
  
	#line_search = cp
    # 线性迭代步的残差下降（相对）量级
 	l_tol = 1e-01
 # l_abs_step_tol = -1e-04
   # 最大线性迭代步	
 	l_max_its = 10
 	
 	# 最大非线性迭代步
 	nl_max_its = 10
 	# 非线性迭代的残值下降（相对）量级
  	#nl_rel_tol = 1e-10
  	# 非线性迭代绝对残值
  	nl_abs_tol = 1e-010

  	
	 abort_on_solve_fail = true	
[]

[Functions]
  [./exact_rho]
    type = CouetteFlowExact
  [../]
[]

[Postprocessors]
  [./l2_err]
    type = ElementL2Error
    variable = rho
    function = exact_rho
  [../]

  [./residuals]
    type = CFDResidual
		execute_on = TIMESTEP
  [../]

  [./elementMaxTimeDerivative]
    type = ElementExtremeTimeDerivative
    variable = rho
	 	execute_on = TIMESTEP_BEGIN
  [../]

  [./area]
    type = AreaPostprocessor
		boundary = left
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
		#verbose = true
    	#setup_log_early = true
    	#time_precision = 6
    	#fit_mode = 100
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
		boundary = 'left right bottom top front back'
		type =NSBC
		variable = rho
	[../]		
	[./x-momentumum_bc]
		boundary = 'left right bottom top front back'
		type =NSBC
		variable = momentum_x
	[../]	
	[./y-momentumum_bc]
		boundary = 'left right bottom top front back'
		type =NSBC
		variable = momentum_y
	[../]
	[./z-momentumum_bc]
		boundary = 'left right bottom top front back'
		type =NSBC
		variable = momentum_z
	[../]		
	[./total-energy_bc]
		boundary = 'left right bottom top front back'
		type =NSBC
		variable = rhoe
	[../]
[]

# 材料属性
[Materials]
  [./cell_materical]
		block = 0
    type = NSCellMaterial
  [../]

  [./face_materical]
		block = 0
    type = NSFaceMaterial
  [../]

  [./bnd_materical]
		boundary = 'left right bottom top front back'
    type = CouetteFlowBndMaterial
  [../]
[]

