[CoupledProblems]
  [./steady]
    input_file = rae2822_wall_distance.i
  [../]
  [./transient]
    input_file = rae2822_sa.i
    [./potential]
      from = steady
      var_name = u
    [../]
  [../]
[]

[Executioner]
  type = SteadyTransientExecutioner
[]
