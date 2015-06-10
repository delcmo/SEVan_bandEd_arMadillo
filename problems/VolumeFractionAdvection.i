# MESH
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmin = 0
  xmax = 1
  block_id = '0'
[]

# GLOBAL PARAMETERS
[GlobalParams]
  # initial Conditions
  pressure_init_left = 1.e6
  pressure_init_right = 1.e6
  vel_init_left = 0.
  vel_init_right = 0.
  temp_init_left = 453.
  temp_init_right = 453.
  liq_vf_init_left = 0.6
  liq_vf_init_right = 0.4
  membrane_position = 0.5
[]

# USEROBJECTS
[UserObjects]
  # liquid
  [./eos_liq]
    type = StiffenedGasEquationOfState
    gamma = 2.35
    p_inf = 1.e9
    q = -1167e3
    cv = 1816.
    q_prime = 0.
  [../]

  # gas phase
  [./eos_gas]
    type = StiffenedGasEquationOfState
    gamma = 1.34
    p_inf = 0
    q = 1968e3
    cv = 1265
  	q_prime = -23e2
  [../]
[]

# VARIABLES
[Variables]
  # liquid phase
  [./alA_liq]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
      type = SbaICs
      eos = eos_liq
    [../]
  [../]

  [./alrhoA_liq]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
      type = SbaICs
      eos = eos_liq
    [../]
  [../]

  [./alrhouA_liq]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
      type = SbaICs
      eos = eos_liq
    [../]
  [../]

  [./alrhoEA_liq]
    family = LAGRANGE
    scaling = 1e-6
    [./InitialCondition]
      type = SbaICs
      eos = eos_liq
    [../]
  [../]

  # gas phase
  [./alrhoA_gas]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
      type = SbaICs
      eos = eos_gas
      isLiquid = false
    [../]
  [../]

  [./alrhouA_gas]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
      type = SbaICs
      eos = eos_gas
      isLiquid = false
    [../]
  [../]

  [./alrhoEA_gas]
    family = LAGRANGE
    scaling = 1e-6
    [./InitialCondition]
      type = SbaICs
      eos = eos_gas
      isLiquid = false
    [../]
  [../]
[]

# KERNELS
[Kernels]
  # liquid volume fraction
  [./VolumedFractionLiqTime]
    type = TimeDerivative
    variable = alA_liq
  [../]

  [./VolumeFractionLiquidConvection]
    type = SbaVolumeFraction
    variable = alA_liq
    alrhoA_k = alrhoA_liq
    alrhouA_x_k = alrhouA_liq
    alrhoEA_k = alrhoEA_liq
    alrhoA_j = alrhoA_gas
    alrhouA_x_j = alrhouA_gas
    alrhoEA_j = alrhoEA_gas
    liquid_volume_fraction = vf_aux_liq
    eos_k = eos_liq
    eos_j = eos_gas
  [../]

  [./VolumeFractionLiquidDissipation]
    type = SbaArtificialDissipation
    variable = alA_liq
    equation_name = VOLUME_FRACTION
    density = density_aux_liq
    pressure = pressure_aux_liq
    velocity_x =   velocity_x_aux_liq
    internal_energy = internal_energy_aux_liq
    liquid_volume_fraction = vf_aux_liq
    eos = eos_liq
  [../]

  # liquid phase (Euler equations)
  [./MassLiquidTime]
    type = TimeDerivative
    variable = alrhoA_liq
  [../]

  [./MomentumLiquidTime]
    type = TimeDerivative
    variable = alrhouA_liq
  [../]

  [./EnergyTimeLiquid]
    type = TimeDerivative
    variable = alrhoEA_liq
  [../]

  # vapor phase
  [./MassTimeGas]
    type = TimeDerivative
    variable = alrhoA_gas
  [../]

  [./MomentumGasTime]
    type = TimeDerivative
    variable = alrhouA_gas
  [../]

  [./EnergyGasTime]
    type = TimeDerivative
    variable = alrhoEA_gas
  [../]
[]

# AUXILARY VARIABLES
[AuxVariables]
  # nodal variables liquid phase
  [./vf_aux_liq]
    family = LAGRANGE
  [../]

  [./velocity_x_aux_liq]
    family = LAGRANGE
  [../]

  [./density_aux_liq]
    family = LAGRANGE
  [../]

  [./internal_energy_aux_liq]
    family = LAGRANGE
  [../]

  [./pressure_aux_liq]
    family = LAGRANGE
  [../]

  # nodal variables gas phase
  [./velocity_x_aux_gas]
    family = LAGRANGE
  [../]

  [./density_aux_gas]
    family = LAGRANGE
  [../]

  [./internal_energy_aux_gas]
    family = LAGRANGE
  [../]

  [./pressure_aux_gas]
    family = LAGRANGE
  [../]

  # elemental variables
  [./beta_max_aux_liq]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./beta_aux_liq]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./beta_max_aux_gas]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./beta_aux_gas]
    family = MONOMIAL
    order = CONSTANT
  [../]

[]

# AUXILARY KERNELS
[AuxKernels]
  # nodal variables liquid phase
  [./VolumeFractionLiquidAK]
    type = VolumeFractionAux
    variable = vf_aux_liq
    alA = alA_liq
  [../]

  [./VelocityLiquidAK]
    type = VelocityAux
    variable = velocity_x_aux_liq
    alrhoA = alrhoA_liq
    alrhouA_x = alrhouA_liq
  [../]

  [./DensitLiquidAK]
    type = DensityAux
    alA = alA_liq
    variable = density_aux_liq
    alrhoA = alrhoA_liq
  [../]

  [./InternalEnergyLiquidAK]
    type = InternalEnergyAux
    variable = internal_energy_aux_liq
    alA = alA_liq
    alrhoA = alrhoA_liq
    alrhouA_x = alrhouA_liq
    alrhoEA = alrhoEA_liq
  [../]

  [./PressureLiquidAK]
    type = PressureAux
    variable = pressure_aux_liq
    alA = alA_liq
    alrhoA = alrhoA_liq
    alrhouA_x = alrhouA_liq
    alrhoEA = alrhoEA_liq
    eos = eos_liq
  [../]

  [./BetaMaxLiquidAK]
    type = MaterialRealAux
    variable = beta_max_aux_liq
    property = beta_max_liq
  [../]

  # elemental variables liquid phase
  [./BetaLiquidAK]
    type = MaterialRealAux
    variable = beta_aux_liq
    property = beta_liq
  [../]

  # Nodal variables gas phase
  [./VelocityGasAK]
    type = VelocityAux
    variable = velocity_x_aux_gas
    alrhoA = alrhoA_gas
    alrhouA_x = alrhouA_gas
  [../]

  [./DensityGasAK]
    type = DensityAux
    variable = density_aux_gas
    alA = alA_liq
    alrhoA = alrhoA_gas
    isLiquid = false
  [../]

  [./InternalEnergyGasAK]
    type = InternalEnergyAux
    variable = internal_energy_aux_gas
    alA = alA_liq
    alrhoA = alrhoA_gas
    alrhouA_x = alrhouA_gas
    alrhoEA = alrhoEA_gas
    isLiquid = false
  [../]

  [./PressureGasdAK]
    type = PressureAux
    variable = pressure_aux_gas
    alA = alA_liq
    alrhoA = alrhoA_gas
    alrhouA_x = alrhouA_gas
    alrhoEA = alrhoEA_gas
    eos = eos_gas
    isLiquid = false
  [../]

  # elemental variables gas phase
  [./BetaMaxGasAK]
    type = MaterialRealAux
    variable = beta_max_aux_gas
    property = beta_max_gas
  [../]

  [./BetaGasAK]
    type = MaterialRealAux
    variable = beta_aux_gas
    property = beta_gas
  [../]
[]


# MATERIALS
[Materials]
  # materials for liquid phase
  [./ViscosityCoefficientsLiquid]
    type = ComputeViscosityCoefficient
    block = '0'
    alrhoA_k = alrhoA_liq
    alrhouA_x_k = alrhouA_liq
    pressure = pressure_aux_liq
    density = density_aux_liq
    volume_fraction_liquid = vf_aux_liq
    eos = eos_liq
  [../]

  # materials for gas phase
  [./ViscosityCoefficientsGas]
    type = ComputeViscosityCoefficient
    block = '0'
    alrhoA_k = alrhoA_gas
    alrhouA_x_k = alrhouA_gas
    pressure = pressure_aux_gas
    density = density_aux_gas
    volume_fraction_liquid = vf_aux_liq
    eos = eos_gas
    isLiquid = false
  [../]

  # materials for interfacial area
  [./InterfacialRelaxationParameter]
    type = InterfacialRelaxationTransfer
    block = '0'
    alrhoA_k = alrhoA_liq
    alrhouA_x_k = alrhouA_liq
    alrhoEA_k = alrhoEA_liq
    alrhoA_j = alrhoA_gas
    alrhouA_x_j = alrhouA_gas
    alrhoEA_j = alrhoEA_gas
    volume_fraction_phase_k = vf_aux_liq
    Aint_max_press = 0.
    Aint_max_vel = 0.
    eos_k = eos_liq
    eos_j = eos_gas
  [../]
[]

# POSTPROCESSORS
[Postprocessors]
[]

# BOUNDARY CONDITIONS
[BCs]
  # liquid phase
  [./VoidFractionLeftLiq]
    type = DirichletBC
    variable = alA_liq
    value = 0.6
    boundary = 'left'
  [../]

  [./MassLiquid]
    type = SbaDirichletBC
    variable = alrhoA_liq
    eos = eos_liq
    boundary = 'left right'
  [../]

  [./MomentumLiquid]
    type = SbaDirichletBC
    variable = alrhouA_liq
    eos = eos_liq
    boundary = 'left right'
  [../]

  [./EnergyLiquid]
    type = SbaDirichletBC
    variable = alrhoEA_liq
    eos = eos_liq
    boundary = 'left right'
  [../]

  # gas phase
  [./MassGas]
    type = SbaDirichletBC
    variable = alrhoA_gas
    eos = eos_gas
    isLiquid = false
    boundary = 'left right'
  [../]

  [./MomentumGas]
    type = SbaDirichletBC
    variable = alrhouA_gas
    eos = eos_gas
    isLiquid = false
    boundary = 'left right'
  [../]

  [./EnergyGas]
    type = SbaDirichletBC
    variable = alrhoEA_gas
    eos = eos_gas
    isLiquid = false
    boundary = 'left right'
  [../]
[]

# PRECONDITIONER
[Preconditioning]
  [./FDP]
    type = FDP
    full = true
    solve_type = 'PJFNK'
    petsc_options_iname = '-mat_fd_coloring_err  -mat_fd_type  -mat_mffd_type'
    petsc_options_value = '1.e-10       ds             ds'
    line_search = 'default'
  [../]
[]

# EXECUTIONER
[Executioner]
  type = Transient
  scheme = 'bdf2'
  num_steps = 2
  end_time = 1.
  dt = 1.e-2
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-7
  l_max_its = 50
  nl_max_its = 10

  [./TimeStepper]
    type = FunctionDT
    time_t =  '0.      1.e-2    2.e-2   4.e-2   0.56'
    time_dt = '1.e-4   1.e-4    1.e-3   1.e-3   1.e-3'
  [../]

  [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]
[]

# OUTPUT
[Outputs]
  output_initial = true
  interval = 1
  console = true
  exodus = true
[]

# DEBUG
[Debug]
  show_var_residual_norms = true
[]
