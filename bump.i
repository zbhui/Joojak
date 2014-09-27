# 全局变量
[GlobalParams]
 	order = FIRST
 	family = MONOMIAL
  	
  mach = 0.5
	reynolds = 10
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
[]

# 网格
[Mesh]
  type = FileMesh
  file = ../high-order-workshop/C1.1_bump/bump_tri_level2.msh
  dim = 2
  
  block_id = 5
  block_name = 'fluid'
  
  boundary_id = '1 2 3 4'
  boundary_name = 'wall top left_in right_out'

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

    petsc_options_iname = '-ksp_type  -pc_type'
  	petsc_options_value = 'gmres 				lu'
	[../]

[]
# 非线性系统求解
[Executioner]
  type = Transient
  solve_type = newton
 	#scheme = 'bdf2'
  num_steps = 1000
  
 	l_tol = 1e-01
 # l_abs_step_tol = -1e-04
 	l_max_its = 50
 	
 	nl_max_its = 20
 	nl_rel_tol = 1e-02
  #nl_abs_tol = 1e-05

	#abort_on_solve_fail = true	
  #end_time = 0.1
  
	[./TimeStepper]
		type = RatioTimeStepper
		dt = 0.100000
		ratio = 2
		step = 2
		max_dt = 1000000
	[../]
[]


[Postprocessors]
	[./area_domain]
		type  = VolumePostprocessor	
	[../]
  [./entropy_err]
    type = BumpElementL2Error
		area = area_domain
  [../]

  [./residuals]
    type = Residual
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
		block = 5
    type = EulerCellMaterial
  [../]

  [./face_material]
		block = 5
    type = EulerFaceMaterial
  [../]

 	[./far_field_material]
		boundary = '2 3 4'
		bc_type = far_field
    type = EulerBndMaterial
  [../]

  [./wall_material]
		boundary = wall
		bc_type = wall
    type = EulerBndMaterial
  [../]

[]

