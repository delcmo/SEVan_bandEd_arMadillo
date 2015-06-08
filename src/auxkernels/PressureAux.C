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
/**
Computes the pressure for a 1-D mesh only.
**/
#include "PressureAux.h"

template<>
InputParameters validParams<PressureAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables:
  params.addRequiredCoupledVar("alA", "liquid alpha*A");  
  params.addRequiredCoupledVar("alrhoA", "alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x", "alpha*rho*u*A_x");
  params.addRequiredCoupledVar("alrhoEA", "alpha*rho*E*A");
  // Aux variables:
  params.addCoupledVar("area", 1., "area");
  // Userobject:
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  // Boolean for phase:
  params.addParam<bool>("isLiquid", true, "is the fluid liquid or not");

  return params;
}

PressureAux::PressureAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variables
    _alA(coupledValue("alA")),
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA_x(coupledValue("alrhouA_x")),
    _alrhoEA(coupledValue("alrhoEA")),
    // Aux variables:
    _area(coupledValue("area")),
    // User Objects for eos
    _eos(getUserObject<EquationOfState>("eos")),
    // Parameter:
    _isLiquid(getParam<bool>("isLiquid"))

{
  mooseAssert(_mesh.dimension() != 1, "The function '"<<_name<<"' can only be used with 1-D mesh.");
}

Real
PressureAux::computeValue()
{
  // Compute phasic volume fraction*A
  Real alA = _isLiquid ? _alA[_qp] : _area[_qp]-_alA[_qp];

  // Compute density, velocity and total energy
  Real rho = _alrhoA[_qp]/alA;
  Real vel = _alrhouA_x[_qp]/_alrhoA[_qp];
  Real rhoE = _alrhoEA[_qp]/alA;

  // Return pressure value
  return _eos.pressure(rho, vel, rhoE);
}
