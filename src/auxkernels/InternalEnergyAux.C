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
This function computes the fluid internal energy 'rhoe' from the conservative variables.
**/
#include "InternalEnergyAux.h"

template<>
InputParameters validParams<InternalEnergyAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("alA", "liquid volume fraction: alpha*A");  
  params.addRequiredCoupledVar("alrhoA", "fluid density: alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x", "fluid x momentum component");
  params.addRequiredCoupledVar("alrhoEA", "alpha*rho*E*A");
  // Coupled aux variables
  params.addCoupledVar("area", 1., "area");
  // Parameters
  params.addParam<bool>("isLiquid", true,"is the flud liquid");

  return params;
}

InternalEnergyAux::InternalEnergyAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
  // Coupled variables:
  _alA(coupledValue("alA")),
  _alrhoA(coupledValue("alrhoA")),
  _alrhouA_x(coupledValue("alrhouA_x")),
  _alrhoEA(coupledValue("alrhoEA")),
  // Aux variables:
  _area(coupledValue("area")),
  // Parameters:
  _isLiquid(getParam<bool>("isLiquid"))
{}

Real
InternalEnergyAux::computeValue()
{
  // Compute the phase volume fraction:
  Real alA = _isLiquid ? _alA[_qp] : _area[_qp]-_alA[_qp];

  // Compute density, velocity and total energy:
  Real rho = _alrhoA[_qp] /alA;
  Real rhoE = _alrhoEA[_qp] / alA;
  Real vel = _alrhouA_x[_qp] / _alrhoA[_qp];

  // Return internal energy:
  return rhoE - 0.5*rho*vel*vel;
}
