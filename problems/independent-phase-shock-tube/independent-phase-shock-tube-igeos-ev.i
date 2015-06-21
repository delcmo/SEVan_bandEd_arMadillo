# MESH
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 200
  xmin = 0
  xmax = 1
  block_id = '0'
[]

# GLOBAL PARAMETERS
[GlobalParams]
  # initial Conditions
  pressure_init_left = 1.e6
  pressure_init_right = 0.1e6
  vel_init_left = 0.
  vel_init_right = 0.
#  rho_init_left = 1.
#  rho_init_right = 0.125
#  temp_init_left = 453.
#  temp_init_right = 453.
  liq_vf_init_left = 0.5
  liq_vf_init_right = 0.5
  membrane_position = 0.5

  # interfacial variables
  interfacial_definition_name = BERRY
  interfacial_variables_on = true

  # stabilization parameters
  is_jump_on = true
  is_first_order_visc = false
  Cjump = 4.
  Ce = 1.
  # cfl
  cfl = 1.
[]

# USEROBJECTS
[UserObjects]
  # eos liquid phase
  [./eos_liq]
    type = StiffenedGasEquationOfState
    gamma = 3.
    p_inf = 0.
    q = 0.
    cv = 0.5
    q_prime = 0.
  [../]

  # eos gas phase
  [./eos_gas]
    type = StiffenedGasEquationOfState
    gamma = 1.4
    p_inf = 0
    q = 0.
    cv = 2.5
  	q_prime = 0.
  [../]

# normalization parmeter liquid
[./norm_param_liq]
  type = NormalizationParameter
  M_threshold = 0.
  a = 0.
  funct_type = Shock_fnct
  velocity_pps_name = AverageVelocityLiq
[../]

  [./norm_param_gas]
    type = NormalizationParameter
    M_threshold = 0.
    a = 0.
    funct_type = Shock_fnct
    velocity_pps_name = AverageVelocityGas
  [../]

# normalization parameter gas

  # jump liquid
  [./jump_press_liq_uo]
    type = JumpGradientInterface
    variable = pressure_aux_liq
    jump_name = jump_press_liq
  [../]

  [./jump_dens_liq_uo]
    type = JumpGradientInterface
    variable = density_aux_liq
    jump_name = jump_dens_liq
  [../]

  # jump gas
  [./jump_press_gas_uo]
    type = JumpGradientInterface
    variable = pressure_aux_gas
    jump_name = jump_press_gas
  [../]

  [./jump_dens_gas_uo]
    type = JumpGradientInterface
    variable = density_aux_gas
    jump_name = jump_dens_gas
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
      rho_init_left = 10.
      rho_init_right = 10.
      eos = eos_liq
    [../]
  [../]

  [./alrhoA_liq]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
      type = SbaICs
      rho_init_left = 10.
      rho_init_right = 10.
      eos = eos_liq
    [../]
  [../]

  [./alrhouA_liq]
    family = LAGRANGE
    scaling = 1e-2
    [./InitialCondition]
      type = SbaICs
      rho_init_left = 10.
      rho_init_right = 10.
      eos = eos_liq
    [../]
  [../]

  [./alrhoEA_liq]
    family = LAGRANGE
    scaling = 1e-4
    [./InitialCondition]
      type = SbaICs
      rho_init_left = 10.
      rho_init_right = 10.
      eos = eos_liq
    [../]
  [../]

  # gas phase
  [./alrhoA_gas]
    family = LAGRANGE
    scaling = 1e+0
    [./InitialCondition]
      type = SbaICs
      rho_init_left = 1.
      rho_init_right = 1.
      eos = eos_gas
      isLiquid = false
    [../]
  [../]

  [./alrhouA_gas]
    family = LAGRANGE
    scaling = 1e-2
    [./InitialCondition]
      type = SbaICs
      rho_init_left = 1.
      rho_init_right = 1.
      eos = eos_gas
      isLiquid = false
    [../]
  [../]

  [./alrhoEA_gas]
    family = LAGRANGE
    scaling = 1e-4
    [./InitialCondition]
      type = SbaICs
      rho_init_left = 1.
      rho_init_right = 1.
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
  [../]

  # liquid phase (Euler equations)
  [./MassLiquidTime]
    type = TimeDerivative
    variable = alrhoA_liq
  [../]
  
  [./MassLiquidAdvection]
    type = SbaMass
    variable = alrhoA_liq
    alrhouA_x_k = alrhouA_liq
  [../]

  [./MassLiquidDissipation]
    type = SbaArtificialDissipation
    variable = alrhoA_liq
    equation_name = CONTINUITY
    density = density_aux_liq
    pressure = pressure_aux_liq
    velocity_x =   velocity_x_aux_liq
    internal_energy = internal_energy_aux_liq
    liquid_volume_fraction = vf_aux_liq
  [../]

  [./MomentumLiquidTime]
    type = TimeDerivative
    variable = alrhouA_liq
  [../]
  
  [./MomentumLiquidAdvection]
    type = SbaMomentum
    variable = alrhouA_liq
    alrhoA_k = alrhoA_liq
    alrhouA_x_k = alrhouA_liq
    alrhoEA_k = alrhoEA_liq
    alrhoA_j = alrhoA_gas
    alrhouA_x_j = alrhouA_gas
    liquid_volume_fraction = vf_aux_liq
    eos_k = eos_liq
  [../]

  [./MomentumLiquidDissipation]
    type = SbaArtificialDissipation
    variable = alrhouA_liq
    equation_name = XMOMENTUM
    density = density_aux_liq
    pressure = pressure_aux_liq
    velocity_x =   velocity_x_aux_liq
    internal_energy = internal_energy_aux_liq
    liquid_volume_fraction = vf_aux_liq
  [../]

  [./EnergyLiquidTime]
    type = TimeDerivative
    variable = alrhoEA_liq
  [../]

  [./EnergyLiquidAdvection]
    type = SbaEnergy
    variable = alrhoEA_liq
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

  [./EnergyLiquidDissipation]
    type = SbaArtificialDissipation
    variable = alrhoEA_liq
    equation_name = ENERGY
    density = density_aux_liq
    pressure = pressure_aux_liq
    velocity_x =   velocity_x_aux_liq
    internal_energy = internal_energy_aux_liq
    liquid_volume_fraction = vf_aux_liq
  [../]

  # vapor phase
  [./MassTimeGas]
    type = TimeDerivative
    variable = alrhoA_gas
  [../]
  
  [./MassGasAdvection]
    type = SbaMass
    variable = alrhoA_gas
    alrhouA_x_k = alrhouA_gas
  [../]

  [./MassGasDissipation]
    type = SbaArtificialDissipation
    variable = alrhoA_gas
    equation_name = CONTINUITY
    density = density_aux_gas
    pressure = pressure_aux_gas
    velocity_x =   velocity_x_aux_gas
    internal_energy = internal_energy_aux_gas
    liquid_volume_fraction = vf_aux_liq
    isLiquid = false
  [../]

  [./MomentumGasTime]
    type = TimeDerivative
    variable = alrhouA_gas
  [../]

  [./MomentumGasAdvection]
    type = SbaMomentum
    variable = alrhouA_gas
    alrhoA_k = alrhoA_gas
    alrhouA_x_k = alrhouA_gas
    alrhoEA_k = alrhoEA_gas
    alrhoA_j = alrhoA_liq
    alrhouA_x_j = alrhouA_liq
    liquid_volume_fraction = vf_aux_liq
    eos_k = eos_gas
    isLiquid = false
  [../]

  [./MomentumGasDissipation]
    type = SbaArtificialDissipation
    variable = alrhouA_gas
    equation_name = XMOMENTUM
    density = density_aux_gas
    pressure = pressure_aux_gas
    velocity_x =   velocity_x_aux_gas
    internal_energy = internal_energy_aux_gas
    liquid_volume_fraction = vf_aux_liq
    isLiquid = false
  [../]

  [./EnergyGasTime]
    type = TimeDerivative
    variable = alrhoEA_gas
  [../]

  [./EnergyGasAdvection]
    type = SbaEnergy
    variable = alrhoEA_gas
    alrhoA_k = alrhoA_gas
    alrhouA_x_k = alrhouA_gas
    alrhoEA_k = alrhoEA_gas
    alrhoA_j = alrhoA_liq
    alrhouA_x_j = alrhouA_liq
    alrhoEA_j = alrhoEA_liq
    liquid_volume_fraction = vf_aux_liq
    eos_k = eos_gas
    eos_j = eos_liq
    isLiquid = false
  [../]

  [./EnergyGasDissipation]
    type = SbaArtificialDissipation
    variable = alrhoEA_gas
    equation_name = ENERGY
    density = density_aux_gas
    pressure = pressure_aux_gas
    velocity_x =   velocity_x_aux_gas
    internal_energy = internal_energy_aux_gas
    liquid_volume_fraction = vf_aux_liq
    isLiquid = false
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

  [./ent_vf_aux_liq]
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

  # elemental variables liquid phase
  [./beta_max_aux_liq]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./beta_aux_liq]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./visc_max_aux_liq]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_aux_liq]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_press_liq]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_dens_liq]
    family = MONOMIAL
    order = CONSTANT
  [../]

  # elemental variables gas phase
  [./beta_aux_gas]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./beta_max_aux_gas]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./visc_max_aux_gas]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./kappa_aux_gas]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_press_gas]
    family = MONOMIAL
    order = CONSTANT
  [../]

  [./jump_dens_gas]
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

  [./EntropyVFLiquidAK]
    type = EntropyVolumeFractionAux
    variable = ent_vf_aux_liq
    volume_fraction_liquid = vf_aux_liq
  [../]

  # elemental variables liquid phase
  [./BetaMaxLiquidAK]
    type = MaterialRealAux
    variable = beta_max_aux_liq
    property = beta_max_liq
  [../]

  [./BetaLiquidAK]
    type = MaterialRealAux
    variable = beta_aux_liq
    property = beta_liq
  [../]
  
  [./KappaMaxLiquidAK]
    type = MaterialRealAux
    variable = visc_max_aux_liq
    property = visc_max_liq
  [../]

  [./KappaLiquidAK]
    type = MaterialRealAux
    variable = kappa_aux_liq
    property = kappa_liq
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
  
  [./KappaMaxGasAK]
    type = MaterialRealAux
    variable = visc_max_aux_gas
    property = visc_max_gas
  [../]

  [./KappaGasAK]
    type = MaterialRealAux
    variable = kappa_aux_gas
    property = kappa_gas
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
    entropy_vf_liquid = ent_vf_aux_liq
    jump_grad_press = jump_press_liq
    jump_grad_rho = jump_dens_liq
    eos = eos_liq
    norm_param = norm_param_liq
    vf_pps_name = InfiniteNormFromAverageValue
  [../]

  # materials for gas phase
  [./ViscosityCoefficientsGas]
    type = ComputeViscosityCoefficient
    block = '0'
    alrhoA_k = alrhoA_gas
    alrhouA_x_k = alrhouA_gas
    pressure = pressure_aux_gas
    density = density_aux_gas
    entropy_vf_liquid = ent_vf_aux_liq
    jump_grad_press = jump_press_gas
    jump_grad_rho = jump_dens_gas
    eos = eos_gas
    norm_param = norm_param_gas
    is_liquid = false
    vf_pps_name = InfiniteNormFromAverageValue
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
    eos_k = eos_liq
    eos_j = eos_gas
  [../]
[]

# POSTPROCESSORS
[Postprocessors]
  [./TimeStepCFL]
    type = TimeStepCFL
    alA_k = alA_liq
    alrhoA_k = alrhoA_liq
    alrhouA_x_k = alrhouA_liq
    alrhoEA_k = alrhoEA_liq
    alrhoA_j = alrhoA_gas
    alrhouA_x_j = alrhouA_gas
    alrhoEA_j = alrhoEA_gas
    eos_k = eos_liq
    eos_j = eos_gas
  [../]

  [./AverageVelocityLiq]
    type = ElementAverageValue
    variable = velocity_x_aux_liq
  [../]

  [./AverageVelocityGas]
    type = ElementAverageValue
    variable = velocity_x_aux_gas
  [../]

  [./AverageVolumeFraction]
    type = ElementAverageValue
    variable = ent_vf_aux_liq
  [../]

  [./InfiniteNormFromAverageValue]
    type = InfiniteNormFromAverageValue
    variable = ent_vf_aux_liq
    name_pps_average = AverageVolumeFraction
  [../]
[]

# BOUNDARY CONDITIONS
[BCs]
  # liquid phase
  [./VoidFractionLeftLiq]
    type = SbaDirichletBC
    variable = alA_liq
    rho_init_left = 10.
    rho_init_right = 10.
    eos = eos_liq
    boundary = 'left'
  [../]

  [./MassLiquidLeft]
    type = SbaDirichletBC
    variable = alrhoA_liq
    rho_init_left = 10.
    rho_init_right = 10.
    eos = eos_liq
    boundary = 'left'
  [../]

  [./MassLiquidRight]
    type = SbaDirichletBC
    variable = alrhoA_liq
    rho_init_left = 10.
    rho_init_right = 10.
    eos = eos_liq
    boundary = 'right'
    leftBoundary = false
  [../]

  [./MomentumLiquidLeft]
    type = SbaDirichletBC
    variable = alrhouA_liq
    rho_init_left = 10.
    rho_init_right = 10.
    eos = eos_liq
    boundary = 'left'
  [../]

  [./MomentumLiquidRight]
    type = SbaDirichletBC
    variable = alrhouA_liq
    rho_init_left = 10.
    rho_init_right = 10.
    eos = eos_liq
    boundary = 'right'
    leftBoundary = false
  [../]

  [./EnergyLiquidLeft]
    type = SbaDirichletBC
    variable = alrhoEA_liq
    rho_init_left = 10.
    rho_init_right = 10.
    eos = eos_liq
    boundary = 'left'
  [../]

  [./EnergyLiquidRight]
    type = SbaDirichletBC
    variable = alrhoEA_liq
    rho_init_left = 10.
    rho_init_right = 10.
    eos = eos_liq
    boundary = 'right'
    leftBoundary = false
  [../]

  # gas phase
  [./MassGasLeft]
    type = SbaDirichletBC
    variable = alrhoA_gas
    rho_init_left = 1.
    rho_init_right = 1.
    eos = eos_gas
    isLiquid = false
    boundary = 'left'
  [../]

  [./MassGasRight]
    type = SbaDirichletBC
    variable = alrhoA_gas
    rho_init_left = 1.
    rho_init_right = 1.
    eos = eos_gas
    isLiquid = false
    boundary = 'right'
    leftBoundary = false
  [../]

  [./MomentumGasLeft]
    type = SbaDirichletBC
    variable = alrhouA_gas
    rho_init_left = 1.
    rho_init_right = 1.
    eos = eos_gas
    isLiquid = false
    boundary = 'left'
  [../]

  [./MomentumGasRight]
    type = SbaDirichletBC
    variable = alrhouA_gas
    rho_init_left = 1.
    rho_init_right = 1.
    eos = eos_gas
    isLiquid = false
    boundary = 'right'
    leftBoundary = false
  [../]

  [./EnergyGasLeft]
    type = SbaDirichletBC
    variable = alrhoEA_gas
    rho_init_left = 1.
    rho_init_right = 1.
    eos = eos_gas
    isLiquid = false
    boundary = 'left'
  [../]

  [./EnergyGasRight]
    type = SbaDirichletBC
    variable = alrhoEA_gas
    rho_init_left = 1.
    rho_init_right = 1.
    eos = eos_gas
    isLiquid = false
    boundary = 'right'
    leftBoundary = false
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
#    line_search = 'default' # does not seem to work well with line_search on
  [../]
[]

# EXECUTIONER
[Executioner]
  type = Transient
  scheme = 'bdf2'
  num_steps = 100
  end_time = 473.e-6
  dt = 1.e-2
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-5
  l_max_its = 50
  nl_max_its = 10

  [./TimeStepper]
    type = PostprocessorDT
    postprocessor = TimeStepCFL
    dt = 1.e-6
  [../]

  [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]
[]

# OUTPUT
[Outputs]
  output_initial = true
  output_final = true
  interval = 1
  console = true
  exodus = true
[]

# DEBUG
[Debug]
#  show_var_residual_norms = true
[]
