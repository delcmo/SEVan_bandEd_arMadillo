#include "StiffenedGasEquationOfState.h"

template<>
InputParameters validParams<StiffenedGasEquationOfState>()
{
  InputParameters params = validParams<EquationOfState>();

  params.addRequiredParam<Real>("gamma", "TODO: describe me");
  params.addRequiredParam<Real>("cv", "Specific heat");
  params.addRequiredParam<Real>("q", "TODO: describe me");
  params.addRequiredParam<Real>("p_inf", "TODO: describe me");

  params.addParam<Real>("q_prime", 0, "TODO: describe me");


  return params;
}


StiffenedGasEquationOfState::StiffenedGasEquationOfState(const std::string & name, InputParameters parameters) :
    EquationOfState(name, parameters),
    _gamma(getParam<Real>("gamma")),
    _cv(getParam<Real>("cv")),
    _q(getParam<Real>("q")),
    _q_prime(getParam<Real>("q_prime")),
    _p_inf(getParam<Real>("p_inf"))
{
  _cp = _cv * _gamma;
}

StiffenedGasEquationOfState::~StiffenedGasEquationOfState()
{
}

Real
StiffenedGasEquationOfState::pressure(Real rho, Real rhou, Real rhoE) const
{
  if (rho == 0.0)
    mooseError("Invalid density of 0.0 detected!");

  return (_gamma - 1) * (rhoE - ((rhou * rhou) / (2 * rho)) - rho * _q) - _gamma * _p_inf;

}

Real
StiffenedGasEquationOfState::temperature(Real rho, Real rhou, Real rhoE) const
{
  if (rho == 0.0)
    mooseError("Invalid density of 0.0 detected!");

  return (1 / _cv) * ((rhoE / rho) - ((rhou * rhou) / (2 * rho * rho)) - _q - (_p_inf / rho));
}

Real
StiffenedGasEquationOfState::c2(Real rho, Real rhou, Real rhoE) const
{
  // NOTE: taken from Marco's code (fish)
  return _gamma * (this->pressure(rho, rhou, rhoE)  + _p_inf) / rho;
}

Real StiffenedGasEquationOfState::c2_from_rho_p(Real rho, Real pressure) const
{
  if (rho == 0.0)
    mooseError("Invalid density of 0.0 detected!");

  return _gamma*(pressure+_p_inf)/rho;
}

Real
StiffenedGasEquationOfState::rho_from_p_T(Real pressure, Real temperature, Real) const
{
  if (((_gamma - 1) * _cv * temperature) == 0.0)
    mooseError("Invalid gamma or cv or temperature detected!");

  return (pressure + _p_inf) / ((_gamma - 1) * _cv * temperature);
}

Real
StiffenedGasEquationOfState::e_from_p_rho(Real pressure, Real rho) const
{
  if ((_gamma - 1) * rho == 0.)
    mooseError("Invalid gamma or density detected!");

  return (pressure + _gamma * _p_inf)/((_gamma - 1) * rho) + _q;
}

Real
StiffenedGasEquationOfState::temp_from_p_rho(Real pressure, Real rho) const
{
  if (rho == 0.0)
    mooseError("Invalid density of 0.0 detected!");

  return (pressure + _p_inf) / ((_gamma - 1) * _cv * rho);
}