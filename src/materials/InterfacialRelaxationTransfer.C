#include "InterfacialRelaxationTransfer.h"
/* This function computes the interfacial Relaxation (PI, velI, PI_bar and velI_bar) and the relaxation parameters (mu and lambda) for the 7 equations model. 
 It also computes the wall heat transfer and friction parameters for each phase.*/
template<>
InputParameters validParams<InterfacialRelaxationTransfer>()
{
  InputParameters params = validParams<Material>();

  params.addParam<std::string>("interfacial_definition_name", "NO_RELAXATION", "Choose definition to compute interfacial variables.");
  params.addParam<Real>("xi_ambrosso", "value for definition of interfacial variables");
  // Conservative variables for phase k:
  params.addRequiredCoupledVar("alrhoA_k", " phase k alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x_k", "x component of phase k alpha*rho*u*A");
  params.addRequiredCoupledVar("alrhoEA_k", "phase k alpha*rho*u*A");
  // Conservative variables for phase j:
  params.addRequiredCoupledVar("alrhoA_j", " phase j alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x_j", "x component of phase j alpha*rho*u*A");
  params.addRequiredCoupledVar("alrhoEA_j", "phase j alpha*rho*u*A");
  // Coupled aux variables:
  params.addCoupledVar("area", 1., "cross-section");
  params.addRequiredCoupledVar("volume_fraction_phase_k", "volume fraction of phase k");
  // Maximum specific interfacial area:
  params.addParam<Real>("Aint_max_press", "Maximum specific interfacial area for pressure relaxation coefficient"); // no need of a default value
  params.addParam<Real>("Aint_max_vel", "Maximum specific interfacial area for velocity relaxation coefficient"); // no need to a default value
  // Equation of states:
  params.addRequiredParam<UserObjectName>("eos_k", "Equation of state for phase k");
  params.addRequiredParam<UserObjectName>("eos_j", "Equation of state for phase j");

  return params;
}

InterfacialRelaxationTransfer::InterfacialRelaxationTransfer(const std::string & name, InputParameters parameters) :
    Material(name, parameters),
    // Definition for interfacial variables:
    _interfacial_param_def("BERRY AMBROSSO LIANG NO_RELAXATION", getParam<std::string>("interfacial_definition_name")),
    _xi(isParamValid("xi_ambrosso") ? getParam<Real>("xi_ambrosso") : 0.),
    // Conservative variables for phase k:
    _alrhoA_k(coupledValue("alrhoA_k")),
    _alrhouA_x_k(coupledValue("alrhouA_x_k")),
    _alrhoEA_k(coupledValue("alrhoEA_k")),
    // Conservative variables for phase j:
    _alrhoA_j(coupledValue("alrhoA_j")),
    _alrhouA_x_j(coupledValue("alrhouA_x_j")),
    _alrhoEA_j(coupledValue("alrhoEA_j")),
    // Coupled aux variables
    _vf_k(coupledValue("volume_fraction_phase_k")),
    _grad_vf_k(coupledGradient("volume_fraction_phase_k")),
    _area(coupledValue("area")),
    // Declare interfacial variables:
    _Aint(declareProperty<Real>("normalized_interfacial_area")),
    _PI(declareProperty<Real>("interfacial_pressure")),
    _velI(declareProperty<Real>("interfacial_velocity")),
    _PI_bar(declareProperty<Real>("average_interfacial_pressure")),
    _velI_bar(declareProperty<Real>("average_interfacial_velocity")),
    // Declare relaxation parameters:
    _press_rel_coeff(declareProperty<Real>("pressure_relaxation")),
    _vel_rel_coeff(declareProperty<Real>("velocity_relaxation")),
    // Maximum specific interfacial area:
    _Aint_max_press(getParam<Real>("Aint_max_press")),
    _Aint_max_vel(getParam<Real>("Aint_max_vel")),
    // Equation of states:
    _eos_k(getUserObject<EquationOfState>("eos_k")),
    _eos_j(getUserObject<EquationOfState>("eos_j"))
{
  // boolean to check for consistent values of the maximum specific interfacial areas
  _has_max_specific_interf_area = false;

  if (!parameters.isParamValid("Aint_max_press") && !parameters.isParamValid("Aint_max_vel"))
    mooseWarning("The maximum specific interfacial area values are not specific in '"<<this->name()<<"': the interfacial and relaxation parameters are set to zero by default.");
  if (!parameters.isParamValid("Aint_max_press") && parameters.isParamValid("Aint_max_vel"))
    mooseError("The maximum specific interfacial area value for the pressure relaxation coefficient is not specified in '"<<this->name()<<"'");
  if (parameters.isParamValid("Aint_max_press") && !parameters.isParamValid("Aint_max_vel"))
    mooseError("The maximum specific interfacial area value for the velocity relaxation coefficient is not specified in '"<<this->name()<<"'");
  if (parameters.isParamValid("Aint_max_press") && parameters.isParamValid("Aint_max_vel"))
    _has_max_specific_interf_area = true;

  // Ambrosso model
  if (!parameters.isParamValid("xi_ambrosso") && _interfacial_param_def==AMBROSSO)
    mooseError("A value for 'xi_ambrosse' was not specified when using the option '"<<_interfacial_param_def<<"' to compute the interfacial variables.");
}

void
InterfacialRelaxationTransfer::computeQpProperties()
{
  // Compute volume fraction of phase j
  Real alpha_j = 1.-_vf_k[_qp];

  // Compute velocity values:
  Real vel_k = _alrhouA_x_k[_qp] / _alrhoA_k[_qp];
  Real vel_j = _alrhouA_x_j[_qp] / _alrhoA_j[_qp];

  // Compute densities:
  Real rho_k = _alrhoA_k[_qp] / (_vf_k[_qp]*_area[_qp]);
  Real rho_j = _alrhoA_j[_qp] / (alpha_j*_area[_qp]);

  // Compute total energies:
  Real rhoE_k = _alrhoEA_k[_qp] / (_vf_k[_qp]*_area[_qp]);
  Real rhoE_j = _alrhoEA_j[_qp] / (alpha_j*_area[_qp]);

  // Compute pressures:
  Real pressure_k = _eos_k.pressure(rho_k, rho_k*vel_k, rhoE_k);
  Real pressure_j = _eos_j.pressure(rho_j, rho_j*vel_j, rhoE_j);

  // Compute the speed of sound for each phase:
  Real c2_k = _eos_k.c2(rho_k, rho_k*vel_k, rhoE_k);
  Real c2_j = _eos_j.c2(rho_j, rho_j*vel_j, rhoE_j);

  // Compute the impedences for each phase:
  Real Z_k = rho_k * std::sqrt(c2_k);
  Real Z_j = rho_j * std::sqrt(c2_j);
  Real sum_Z = Z_k + Z_j;

  // Compute sign of volume fraction gradient of phase k
  Real sgn_k = _grad_vf_k[_qp](0) > 0 ? 1. : -1.;

  /*******************************************************************************/
  /*********************** Compute interfacial variables *************************/
  /*******************************************************************************/
  // declare temporary variables used in the definition of the interfacial variables
  Real beta(0.), mu(0.), temp_k(0.), temp_j(0.);

  // Compute the interfacial variables
  switch (_interfacial_param_def)
  {
    case BERRY:
      // Compute the average interfacial Relaxation parameters:
      _PI_bar[_qp] = ( Z_j*pressure_k + Z_k*pressure_j ) / sum_Z;
      _velI_bar[_qp] = ( Z_k*vel_k + Z_j*vel_j ) / sum_Z;
      
      // Compute interfacial Relaxation parameters:
      _PI[_qp] = _PI_bar[_qp] + Z_k*Z_j/sum_Z * (vel_j-vel_k);
      _velI[_qp] = _velI_bar[_qp] + sgn_k * (pressure_j-pressure_k)/sum_Z;
      break;
    case AMBROSSO:
      // Compute interfacial velocity:
      beta = _xi*_vf_k[_qp]*rho_k;
      beta *= 1./(_xi*_vf_k[_qp]*rho_k + (1.-_xi)*(1.-_vf_k[_qp])*rho_j);
      _velI_bar[_qp] = beta*vel_k + (1.-beta)*vel_j;
      _velI[_qp] = _velI_bar[_qp];
      // Compute interfacial pressure:
      temp_k = _eos_k.temp_from_p_rho(pressure_k, rho_k);
      temp_j = _eos_j.temp_from_p_rho(pressure_j, rho_j);
      mu = (1.-beta)*temp_j/(beta*temp_k+(1.-beta)*temp_j);
      _PI_bar[_qp] = mu*pressure_k + (1.-mu)*pressure_j;
      _PI[_qp] = _PI_bar[_qp];
      break;
    case LIANG:
      // Compute intefacial velocity:
      _velI_bar[_qp] = _vf_k[_qp]*rho_k*vel_k + (1.-_vf_k[_qp])*rho_j*vel_j;
      _velI_bar[_qp] *= 1./(_vf_k[_qp]*rho_k + (1.-_vf_k[_qp])*rho_j);
      _velI[_qp] = _velI_bar[_qp];
      // Compute interfacial pressure:
      _PI_bar[_qp] = _vf_k[_qp]*pressure_k + (1.-_vf_k[_qp])*pressure_j;
      _PI[_qp] = _PI_bar[_qp];
      break;
    case NO_RELAXATION:
      // Compute intefacial velocity:
      _velI_bar[_qp] = 0.; _velI[_qp] = 0.;
      // Compute interfacial pressure:
      _PI_bar[_qp] = 0.; _PI[_qp] = 0.;
      break;
    default:
      mooseError("The intercial definition '"<<_interfacial_param_def<<"' is not implemented in the function '"<<this->name()<<"'.");
      break;
  }

  /*******************************************************************************/
  /*********************** Compute relaxation parameters *************************/
  /*******************************************************************************/
  if (_has_max_specific_interf_area)
  {
    _Aint[_qp] = 6.75*(1-_vf_k[_qp])*(1-_vf_k[_qp])*_vf_k[_qp];
    _press_rel_coeff[_qp] = _Aint_max_press*_Aint[_qp]/sum_Z; /*(mu)*/
    _vel_rel_coeff[_qp] = _Aint_max_vel*_Aint[_qp]*Z_j*Z_k/sum_Z; /*(lambda)*/
  }
  else
  {
    _Aint[_qp] = 0.;
    _press_rel_coeff[_qp] = 0.; _vel_rel_coeff[_qp] = 0.;
  }
}