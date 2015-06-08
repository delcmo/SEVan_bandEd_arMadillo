# GLOBAL PARAMETERS
[GlobalParams]
# Other parameters
viscosity_name = FIRST_ORDER

# Mass and heat transfer
isJumpOn = false
isMassOn = false
isHeatOn = false
isShock = false
Aint = 0.

# Initial Conditions
pressure_init_left = 1e6
pressure_init_right = 0.5e6
vel_init_left = 10
vel_init_right = 10
temp_init_left = 453.
temp_init_right = 453.
alpha_init_left = 0.5
alpha_init_right = 0.5
membrane = 0.5
length = 1.
[]

# USER OBJECTS
[UserObjects]
  # liquid
  [./eos_liq]
    type = EquationOfState
    gamma = 2.35
    Pinf = 1.e9
    q = -1167e3
    Cv = 1816.
    q_prime = 0.
  [../]

  # gas phase
  [./eos_gas]
    type = EquationOfState
    gamma = 1.34 # 1.43
    Pinf = 0
    q = 1968e3 # 2030e3
    Cv = 1265 # 1040
  	q_prime = -23e2
  [../]
[]

# MESH
[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = 0
  xmax = 1
  block_id = '0'
[]

# VARIABLES
[Variables]
  # liquid phase
  [./alA_liq]
    family = LAGRANGE
    scaling = 1e-2
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_liq
    [../]
  [../]

  [./alrhoA_liq]
    family = LAGRANGE
    scaling = 1e-2
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_liq
    [../]
  [../]

  [./alrhouA_liq]
    family = LAGRANGE
    scaling = 1e-4
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_liq
    [../]
  [../]

  [./alrhoEA_liq]
    family = LAGRANGE
    scaling = 1e-7
    [./InitialCondition]
          type = ConservativeVariables1DXIC
          area = area
          eos = eos_liq
    [../]
  [../]

  # gas phase
  [./alrhoA_gas]
    family = LAGRANGE
    scaling = 1e-2
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_gas
        isLiquid = false
    [../]
  [../]

  [./alrhouA_gas]
    family = LAGRANGE
    scaling = 1e-4
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
        eos = eos_gas
        isLiquid = false
    [../]
  [../]

  [./alrhoEA_gas]
    family = LAGRANGE
    scaling = 1e-7
    [./InitialCondition]
        type = ConservativeVariables1DXIC
        area = area
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
    alrhouA_k = alrhouA_liq
    alrhoEA_k = alrhoEA_liq
    liquid_volume_fraction = vf_liq
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
    liquid_volume_fraction = vf_liq
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

#  [./MomConvLiq]
#    type = EelMomentum
#    variable = alrhouA_l
#    vel_x = velocity_x_aux_l
#    vel_x_2 = velocity_x_aux_g
#    pressure = pressure_aux_l
#    area = area_aux
#    vf_liquid = alpha_aux_l
#  [../]

#  [./EnergyConvLiq]
#    type = EelEnergy
#    variable = alrhoEA_l
#    alrhoA = alrhoA_l
#    alrhouA_x = alrhouA_l
#    vel_x_2 = velocity_x_aux_g
#    pressure_liq = pressure_aux_l
#    pressure_gas = pressure_aux_g
#    area = area_aux
#    vf_liquid = alpha_aux_l
#    eos = eos_liq
#  [../]

#  [./MassViscLiq]
#    type = EelArtificialVisc
#    variable = alrhoA_l
#    equation_name = CONTINUITY
#    density = density_aux_l
#    pressure = pressure_aux_l
#    velocity_x = velocity_x_aux_l
#    internal_energy = internal_energy_aux_l
#    area = area_aux
#    vf_liquid = alpha_aux_l
#    eos = eos_liq
#  [../]

#  [./MomViscLiq]
#    type = EelArtificialVisc
#    variable = alrhouA_l
#    equation_name = XMOMENTUM
#    density = density_aux_l
#    pressure = pressure_aux_l
#    velocity_x = velocity_x_aux_l
#    internal_energy = internal_energy_aux_l
#    area = area_aux
#    vf_liquid = alpha_aux_l
#    eos  =eos_liq
#  [../]

#  [./EnergyViscLiq]
#    type = EelArtificialVisc
#    variable = alrhoEA_l
#    equation_name = ENERGY
#    density = density_aux_l
#    pressure = pressure_aux_l
#    velocity_x = velocity_x_aux_l
#    internal_energy = internal_energy_aux_l
#    area = area_aux
#    vf_liquid = alpha_aux_l
#    eos = eos_liq
#  [../]

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

#  [./MassConvGas]
#    type = EelMass
#    variable = alrhoA_g
#    alrhouA_x = alrhouA_g
#    area = area_aux
#    isLiquid = false
#  [../]

#  [./MomConvGas]
#    type = EelMomentum
#    variable = alrhouA_g
#    vel_x = velocity_x_aux_g
#    vel_x_2 = velocity_x_aux_l
#    pressure = pressure_aux_g
#    area = area_aux
#    vf_liquid = alpha_aux_l
#    isLiquid = false
#  [../]

#  [./EnergyConvGas]
#    type = EelEnergy
#    variable = alrhoEA_g
#    alrhoA = alrhoA_g
#    alrhouA_x = alrhouA_g
#    vel_x_2 = velocity_x_aux_l
#    pressure_liq = pressure_aux_l
#    pressure_gas = pressure_aux_g
#    area = area_aux
#    vf_liquid = alpha_aux_l
#    isLiquid = false
#    eos = eos_gas
#  [../]

#  [./MassViscGas]
#    type = EelArtificialVisc
#   variable = alrhoA_g
#    equation_name = CONTINUITY
#    density = density_aux_g
#    pressure = pressure_aux_g
#    velocity_x = velocity_x_aux_g
#    internal_energy = internal_energy_aux_g
#    area = area_aux
#    vf_liquid = alpha_aux_l
#    isLiquid = false
#    eos  =eos_gas
#  [../]

#  [./MomViscGas]
#    type = EelArtificialVisc
#    variable = alrhouA_g
#    equation_name = XMOMENTUM
#    density = density_aux_g
#    pressure = pressure_aux_g
#    velocity_x = velocity_x_aux_g
#    internal_energy = internal_energy_aux_g
#    area = area_aux
#    vf_liquid = alpha_aux_l
#    isLiquid = false
#    eos = eos_gas
#  [../]

#  [./EnergyViscGas]
#    type = EelArtificialVisc
#    variable = alrhoEA_g
#    equation_name = ENERGY
#    density = density_aux_g
#    pressure = pressure_aux_g
#    velocity_x = velocity_x_aux_g
#    internal_energy = internal_energy_aux_g
#    area = area_aux
#    vf_liquid = alpha_aux_l
#    isLiquid = false
#    eos = eos_gas
# [../]
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

  # nodel variables gas phase
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
    variable = alpha_aux_liq
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
    variable = density_aux_liq
    alrhoA = alrhoA_liq
    vf_liquid = alpha_aux_liq
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
    alrhouA = alrhouA_gas
  [../]

  [./DensityGasAK]
    type = DensityAux
    variable = density_aux_gas
    alrhoA = alrhoA_gas
    vf_liquid = alpha_aux_l
    area = area_aux
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

  [./PressureLiquidAK]
    type = PressureAux
    variable = pressure_aux_gas
    alA = alA_liq
    alrhoA = alrhoA_gas
    alrhouA_x = alrhouA_gas
    alrhoEA = alrhoEA_gas
    eos = eos_gas
    isLiquid = false
  [../]

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
    alrhouA_k = alrhouA_liq
    pressure = pressure_aux_liq
    density = density_aux_liq
    volume_fraction_liquid = alpha_aux_liq
    eos = eos_liq
  [../]

  # materials for gas phase
  [./ViscosityCoefficientsGas]
    type = ComputeViscosityCoefficient
    block = '0'
    alrhoA_k = alrhoA_gas
    alrhouA_k = alrhouA_gas
    pressure = pressure_aux_gas
    density = density_aux_gas
    volume_fraction_liquid = alpha_aux_liq
    eos = eos_gas
  [../]

  # materials for interfacial area
  [./InterfacialRelaxationParameter]
    type = InterfacialRelaxationTransfer
    block = '0'
    alrhoA_k = alrhoA_liq
    alrhouA_k = alrhouA_liq
    alrhoEA_k = alrhoEA_liq
    alrhoA_j = alrhoA_gas
    alrhouA_j = alrhouA_gas
    alrhoEA_j = alrhoEA_gas
    volume_fraction_phase_k = vf_aux_liq
    eos_liq = eos_liq
    eos_gas = eos_gas
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
    value = 0.5
    boundary = 'left'
  [../]

  [./MassLiquid]
    type = DirichletBC
    variable = alrhoA_liq
    value = 1.
    boundary = 'left right'
  [../]

  [./MomentumLiquid]
    type = DirichletBC
    variable = alrhouA_liq
    value = 0.
    boundary = 'left right'
  [../]

  [./EnergyLiquid]
    type = DirichletBC
    variable = alrhoEA_liq
    value = 1.
    boundary = 'left right'
  [../]

  # gas phase
  [./MassGas]
    type = DirichletBC
    variable = alrhoA_gas
    value = 1.
    boundary = 'left'
  [../]

  [./MomentumGas]
    type = DirichletBC
    variable = alrhouA_gas
    value = 0
    boundary = 'left right'
  [../]

  [./EnergyGas]
    type = DirichletBC
    variable = alrhoEA_gas
    value = 1.
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
#  num_steps = 1
  end_time = 1.
  dt = 1.e-2
  dtmin = 1e-9
  l_tol = 1e-8
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-9
  l_max_its = 50
  nl_max_its = 10
[./TimeStepper]
    type = FunctionDT
#    time_t =  '0.      2.e-4    1.e-2   2.e-2   0.56'
#    time_dt = '1.e-5   1.e-4    1.e-4   1.e-3   1.e-3'
    time_t =  '0.      1.e-2    2.e-2   4.e-2   0.56'
    time_dt = '1.e-4   1.e-4    1.e-3   1.e-3   1.e-3'
#    time_t =  '0.      1.e-2    3.e-2   1.e-1   0.56'
#    time_dt = '1.e-2   1.e-3    1.e-3   1.e-3   1.e-3'
  [../]

  [./Quadrature]
    type = GAUSS
    order = SECOND
  [../]
[]

# OUTPUT
[Outputs]
  output_initial = true
  interval = 20
  console = true
  exodus = true
  postprocessor_screen = true
  perf_log = true
[]

##############################################################################################
#                                        DEBUG                                               #
##############################################################################################
# Debug                 #
##############################################################################################

#[Debug]
#  show_var_residual_norms = true
#[]
