[GlobalParams]
  use_displaced_mesh = true
[]

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
  displacements = 'disp_x disp_y'
[]

[Kernels]
  [./disp_x_time]
    type = EmptyTimeDerivative
    variable = disp_x
  [../]
  [./disp_y_time]
    type = EmptyTimeDerivative
    variable = disp_y
  [../]

  [./disp_x_space]
    type = ElasticityKernel
    variable = disp_x
  [../]
  [./disp_y_space]
    type = ElasticityKernel
    variable = disp_y
  [../]
[]

[BCs]
  [./disp_x_bc]
    type = DirichletBC
    value = 1
    variable = disp_x
    boundary = left
  [../]
  [./disp_y_bc]
    type = DirichletBC
    value = 0
    variable = disp_y
    boundary = left
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  scheme = 'crank-nicolson'
  dt = 0.02
  end_time = 1
  num_steps = 100
  l_tol = 1e-04
  l_max_its = 10
  nl_max_its = 10
  nl_rel_tol = 1e-08
[]


[Outputs]
  [./exodus]
    type = Exodus
  [../]
  [./console]
    type = Console	
    perf_log = true
    output_on = linear
  [../]
[]

[Materials]
  [./cell_materical]
    block = 0
    type = LinearElasticityMaterial
    youngs_modulus = 1
    poissons_ratio = 0.25
    disp_x = disp_x
    disp_y = disp_y
  [../]
[]
