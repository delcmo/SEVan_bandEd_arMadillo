/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "SbaICs.h"

template<>
InputParameters validParams<SbaICs>()
{
  InputParameters params = validParams<InitialCondition>();

  // Name of the function to use to compute the area
  params.addParam<FunctionName>("area", "function to compute the area");
  // Initial pressure conditions:
  params.addRequiredParam<Real>("pressure_init_left", "Left initial value of the pressure");
  params.addRequiredParam<Real>("pressure_init_right", "Right initial value of the pressure");
  // Initial velocity conditions:
  params.addRequiredParam<Real>("vel_init_left", "Left initial value of the velocity");
  params.addRequiredParam<Real>("vel_init_right", "Right Inital value of the velocity");
  // Initial temperature conditions:
  params.addParam<Real>("temp_init_left", "Left initil value of the temperature");
  params.addParam<Real>("temp_init_right", "Right initil value of the temperature");
  // Initial density conditions:
  params.addParam<Real>("rho_init_left_liq", "Left initial value of the density");
  params.addParam<Real>("rho_init_right_liq", "Right initial value of the density");
  // Initial volume fraction conditions
  params.addRequiredParam<Real>("liq_vf_init_left", "Left initial value of the LIQUID volume fraction");
  params.addRequiredParam<Real>("liq_vf_init_right", "Right initial value of the LIQUID volume fraction");
  // Membrane position:
  params.addParam<Real>("membrane_position", 0.5, "The value of the membrane");
  params.addParam<Real>("length", 0.01, "To smooth the IC over a given length");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "parameters for eos.");
  // Boolean
  params.addParam<bool>("isLiquid", true, "is phase liquid or not?");

  return params;
}

SbaICs::SbaICs(const std::string & name,
                     InputParameters parameters) :
    InitialCondition(name, parameters),
    // Cross section function
    _area(getFunction("area")),
    // Initial pressure values:
    _p_left(getParam<Real>("pressure_init_left")),
    _p_right(getParam<Real>("pressure_init_right")),
    // Initial velocity values:
    _v_left(getParam<Real>("vel_init_left")),
    _v_right(getParam<Real>("vel_init_right")),
    // Initial temperature values:
    _t_left(getParam<Real>("temp_init_left")),
    _t_right(getParam<Real>("temp_init_right")),
    // Initial density values:
    _rho_left(getParam<Real>("rho_init_left_liq")),
    _rho_right(getParam<Real>("rho_init_right_liq")),
    // Initial liquid volume fraction values:
    _liq_vf_left(getParam<Real>("alpha_init_left")),
    _liq_vf_right(getParam<Real>("alpha_init_right")),
    // Position of the membrane:
    _membrane(getParam<Real>("membrane_position")),
    _length(getParam<Real>("length")),
    // Equation of State:
    _eos(getUserObject<EquationOfState>("eos")),
    // Boolean:
    _isLiquid(getParam<bool>("isLiquid"))
{
  // Determine the initial condition type:
  if (parameters.isParamValid("rho_init_left") && parameters.isParamValid("vel_init_left") && parameters.isParamValid("press_init_left"))
    _ics_type = 0;
  else if (parameters.isParamValid("temp_init_left") && parameters.isParamValid("vel_init_left") && parameters.isParamValid("press_init_left"))
    _ics_type = 1;
  else if (parameters.isParamValid("rho_init_left") && parameters.isParamValid("vel_init_left") && parameters.isParamValid("temp_init_left"))
    _ics_type = 2;
  else
    mooseError("The input values provided in the input files are incomplete.");

  // Boolean for area function
  _isArea = parameters.isParamValid("area") ? true : false;
}

Real
SbaICs::value(const Point & p)
{
  // Compute the x1 and x2
  Real x1 = _membrane - 0.5 * _length;
  Real x2 = x1 + _length;

  // Compute the density, the momentum and the total energy from the input values:
  Real a_linear, b_linear, vf, rho, rhou, rhoE;

  if (_ics_type==0)
  {
    if (p(0)<x1)
    {
      vf = _isLiquid ? _liq_vf_left : 1.-_liq_vf_left;
      rho = _rho_left;
      rhou = rho*_v_left;
      rhoE = rho*( _eos.e_from_p_rho(_p_left, rho) + 0.5*_v_left*_v_left);
    }
    else
    {
      vf = _isLiquid ? _liq_vf_right : 1.-_liq_vf_right;
      rho = _rho_right;
      rhou = rho*_v_right;
      rhoE = rho*( _eos.e_from_p_rho(_p_right, rho) + 0.5*_v_right*_v_right);
    }
  }
  else if (_ics_type==1)
  {
    if (p(0)<x1)
    {
      vf = _isLiquid ? _liq_vf_left : 1.-_liq_vf_left;
      rho = _eos.rho_from_p_T(_p_left, _t_left);
      rhou = rho*_v_left;
      rhoE = rho*( _eos.e_from_p_rho(_p_left, rho) + 0.5*_v_left*_v_left);
    }
    else
    {
      vf = _isLiquid ? _liq_vf_right : 1.-_liq_vf_right;
      rho = _eos.rho_from_p_T(_p_right, _t_right);
      rhou = rho*_v_right;
      rhoE = rho*( _eos.e_from_p_rho(_p_right, rho) + 0.5*_v_right*_v_right);
    }
  }
  else // _ics_type==2
  {
    if (p(0)<x1)
    {
      vf = _isLiquid ? _liq_vf_left : 1.-_liq_vf_left;
      rho = _rho_left;
      rhou = rho*_v_left;
      Real p = _eos.p_from_rho_T(rho, _t_left);
      rhoE = rho*( _eos.e_from_p_rho(p, rho) + 0.5*_v_left*_v_left);
    }
    else
    {
      vf = _isLiquid ? _liq_vf_right : 1.-_liq_vf_right;
      rho = _rho_right;
      rhou = rho*_v_right;
      Real p = _eos.p_from_rho_T(rho, _t_right);
      rhoE = rho*( _eos.e_from_p_rho(p, rho) + 0.5*_v_right*_v_right);
    }
  }

  // Value of the area:
  Real area = _isArea ? _area.value(0., p) : 1.;

  // String name for phase:
  std::string phase_name = _isLiquid ? "_liq" : "_gas";

  // Return the values
  if ( _var.name() == "alA"+phase_name) // volume fraction
      return vf*area;
  else if ( _var.name() == "alrhoA"+phase_name) // density
    return vf*rho*area;
  else if ( _var.name() == "alrhouA"+phase_name) // momentum
    return vf*rhou*area;
  else if ( _var.name() == "alrhoEA"+phase_name) // energy
    return vf*rhoE*area;
  else
    mooseError("In '"<<this->name()<<"', the variable name '"<<_var.name()<<" cannot be used in this function. The variable names must be chosen among: alA_liq, alrhoA_liq, alrhouA_liq, alrhoEA_liq, alrhoA_gas, alrhouA_gas and alrhoEA_gas");
}