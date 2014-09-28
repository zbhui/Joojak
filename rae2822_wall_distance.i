[Mesh]
  type = FileMesh
  file = ../high-order-workshop/C2.2_rae2822/rae2822_level5.msh
  dim = 2

  boundary_id = '1 3' 
  boundary_name = 'wall far_field'

  block_id = '0'
  block_name = 'fluid'
[]

[BCs]
  [./wall]
    type = DirichletBC
    variable = u
    boundary = wall
    value = 0
  [../]

  [./other]
		type = NeumannBC
		variable = u
		boundary  = 'far_field'
		value = 0
  [../]
[]

[Functions]
  [./forcing_fn]
    type = ParsedFunction
    value = 1
  [../]
[]

[Variables]
  [./u]
		order = FIRST
    family = LAGRANGE
		initial_condition = 0
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
  [./forcing_fn]
    type = UserForcingFunction
    variable = u
    function = forcing_fn
  [../]
[]

[Executioner]
  type = Steady
  l_tol = 1e-4
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-10
  solve_type = NEWTON
[]

[Outputs]
	exodus = true
  console = true
	linear_residuals = true
	nonlinear_residuals = true
[]
