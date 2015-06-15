#include "TimeStepCFL.h"

template<>
InputParameters validParams<TimeStepCFL>()
{
  InputParameters params = validParams<ElementPostprocessor>();

  // Conservative variables for phase k:
  params.addRequiredCoupledVar("alA_k", " phase k alpha*A");
  params.addRequiredCoupledVar("alrhoA_k", " phase k alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x_k", "x component of phase k alpha*rho*u*A");
  params.addRequiredCoupledVar("alrhoEA_k", "phase k alpha*rho*u*A");
  // Conservative variables for phase j:
  params.addRequiredCoupledVar("alrhoA_j", " phase k alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x_j", "x component of phase j alpha*rho*u*A");
  params.addRequiredCoupledVar("alrhoEA_j", "phase j alpha*rho*u*A");
  // Coupled aux variables:
  params.addCoupledVar("area", 1., "cross-section");
  // Equation of states:
  params.addRequiredParam<UserObjectName>("eos_k", "Equation of state for phase k");
  params.addRequiredParam<UserObjectName>("eos_j", "Equation of state for phase j");
  // Parameter
  params.addParam<Real>("cfl", 0.8, "CFL number to supply by the user");

  return params;
}

TimeStepCFL::TimeStepCFL(const std::string & name, InputParameters parameters) :
    ElementPostprocessor(name, parameters),
    // Conservative variables for phase k:
    _alA_k(coupledValue("alA_k")),
    _alrhoA_k(coupledValue("alrhoA_k")),
    _alrhouA_x_k(coupledValue("alrhouA_x_k")),
    _alrhoEA_k(coupledValue("alrhoEA_k")),
    // Conservative variables for phase j:
    _alrhoA_j(coupledValue("alrhoA_j")),
    _alrhouA_x_j(coupledValue("alrhouA_x_j")),
    _alrhoEA_j(coupledValue("alrhoEA_j")),
    // Coupled aux variables
    _area(coupledValue("area")),
    // Equation of states:
    _eos_k(getUserObject<EquationOfState>("eos_k")),
    _eos_j(getUserObject<EquationOfState>("eos_j")),
    // Interfacial velocity
    _velI(getMaterialProperty<Real>("interfacial_velocity")),
    // Parameters
    _cfl(getParam<Real>("cfl")),
    _value(0.)
{
}

TimeStepCFL::~TimeStepCFL()
{
}

void
TimeStepCFL::initialize()
{
  _value = std::numeric_limits<Real>::max();
}

void
TimeStepCFL::execute()
{
  // Compute cell size
  Real h_cell = std::pow(_current_elem->volume(), 1./_mesh.dimension());

  // Loop over quadrature points
  for (unsigned qp = 0; qp < _qrule->n_points(); ++qp)
  {
    // Compute time step for phase k
    Real rho_k = _alrhoA_k[qp] / _alA_k[qp];
    Real vel_k = _alrhouA_x_k[qp] / _alrhoA_k[qp];
    Real rhoE_k = _alrhoEA_k[qp]/_alA_k[qp];
    Real eigen_k = std::fabs(vel_k)+std::sqrt(_eos_k.c2(rho_k, vel_k*rho_k, rhoE_k));
    Real dt_k = _cfl * h_cell / eigen_k;

    // Compute time step for phase j
    Real rho_j = _alrhoA_j[qp] / (_area[qp]-_alA_k[qp]);
    Real vel_j = _alrhouA_x_j[qp] / _alrhoA_j[qp];
    Real rhoE_j = _alrhoEA_j[qp]/(_area[qp]-_alA_k[qp]);
    Real eigen_j = std::fabs(vel_j)+std::sqrt(_eos_j.c2(rho_j, vel_j*rho_j, rhoE_j));
    Real dt_j = _cfl * h_cell / eigen_j;

    // Compute time step for interfacial velocity
    Real dt_int = _cfl * h_cell / std::fabs(_velI[qp]);

    // Compute the local time step
    _value = std::min(std::min(_value, dt_int), std::min(dt_k, dt_j));
  }
}

Real
TimeStepCFL::getValue()
{
  _communicator.min(_value);
  return _value;
}

void
TimeStepCFL::threadJoin(const UserObject & uo)
{
  const TimeStepCFL & pps = dynamic_cast<const TimeStepCFL &>(uo);
  _value = std::min(_value, pps._value);
}