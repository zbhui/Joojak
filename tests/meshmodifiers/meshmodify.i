[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 10
  ny = 10
  nz = 10
[]

[MeshModifiers]
  [./extrude]
    type = BuildSideSetFromBlock
    new_boundary = boundary_from_block
  [../]
[]

[Variables]
  [./u]
  [../]
[]

[AuxVariables]
  [./distance_to_left_nodes]
  [../]
  [./penetration_to_left]
  [../]
[]

[Kernels]
  [./diff]
    type = Diffusion
    variable = u
  [../]
[]

[AuxKernels]
  [./nodal_distance_aux]
    type = NearestNodeDistanceAux
    variable = distance_to_left_nodes
    boundary = boundary_from_block
    paired_boundary = 0
  [../]
  [./penetration_aux]
    type = PenetrationAux
    variable = penetration_to_left
    boundary = boundary_from_block
    paired_boundary = 1
  [../]
[]


[Executioner]
  # Preconditioned JFNK (default)
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  output_initial = true
  exodus = true
  print_linear_residuals = true
  print_perf_log = true
[]
