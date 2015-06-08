#include "EquationOfState.h"
#include "MooseError.h"

template<>
InputParameters validParams<EquationOfState>()
{
  InputParameters params = validParams<UserObject>();

  params.addPrivateParam<MultiMooseEnum>("execute_on");
  params.addPrivateParam<bool>("use_displaced_mesh");
  params.registerBase("EquationOfState");

  return params;
}

EquationOfState::EquationOfState(const std::string & name, InputParameters parameters) :
    GeneralUserObject(name, parameters)
{}

EquationOfState::~EquationOfState()
{
  // Destructor, empty
}

Real EquationOfState::pressure(Real, Real, Real) const
{
  this->error_not_implemented("pressure");
  return 0.;
}

Real EquationOfState::temperature(Real, Real, Real) const
{
  this->error_not_implemented("temperature");
  return 0.;
}

Real
EquationOfState::c2(Real, Real, Real) const
{
  this->error_not_implemented("c2");
  return 0.;
}

Real
EquationOfState::c2_from_rho_p(Real, Real) const
{
  this->error_not_implemented("c2_from_rho_p");
  return 0.;
}

Real EquationOfState::rho_from_p_T(Real, Real) const
{
  this->error_not_implemented("rho_from_p_T");
  return 0.;
}

Real EquationOfState::p_from_rho_T(Real, Real) const
{
  this->error_not_implemented("p_from_rho_T");
  return 0.;
}

Real
EquationOfState::e_from_p_rho(Real, Real) const
{
  this->error_not_implemented("e_from_p_rho");
  return 0.;
}

Real
EquationOfState::temp_from_p_rho(Real, Real) const
{
  this->error_not_implemented("temperature");
  return 0.;
}

void EquationOfState::error_not_implemented(std::string method_name) const
{
  mooseError("Your EquationOfState object does not implement: " + method_name);
}

