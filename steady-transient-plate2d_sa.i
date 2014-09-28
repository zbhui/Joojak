[CoupledProblems]
  [./steady]
    input_file = plate2d_wall_distance.i
  [../]
  [./transient]
    input_file = plate2d_sa.i
    [./potential]
      from = steady
      var_name = u
    [../]
  [../]
[]

[Executioner]
  type = SteadyTransientExecutioner
[]
