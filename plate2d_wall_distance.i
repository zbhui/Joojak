[Mesh]
  type = FileMesh
  file = ../high-order-workshop/C1.4_plate/a2-125-2s.msh
  dim = 2

  boundary_id = '1 2 3 4 5' 
  boundary_name = 'symmetric wall right top left'

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
		boundary  = 'symmetric right top left'
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
  console = true
	linear_residuals = true
	nonlinear_residuals = true
[]
