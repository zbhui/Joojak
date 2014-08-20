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


# 非线性系统求解
[Executioner]
  type = Transient
  solve_type = PJFNK
 	scheme = 'bdf2'
  dt = 0.005
  num_steps = 1
  
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
  
  #petsc_options_iname = '-pc_type -pc_hypre_type'
  #petsc_options_value = 'hypre boomeramg'
  	
	abort_on_solve_fail = true	
  #end_time = 0.1
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

# 输出和后处理
[Outputs]
  file_base = isovortex
	[./exodus]
		type = Exodus
		output_initial = true
	[../]
	
	[./console]
		type = Console	
		perf_log = true
		linear_residuals = true
	  nonlinear_residuals =  true	
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
		type =TimeDerivative
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
