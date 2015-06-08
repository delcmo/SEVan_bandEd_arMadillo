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
This function computes the density of the fluid.
**/
#include "DensityAux.h"

template<>
InputParameters validParams<DensityAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variables
  params.addRequiredCoupledVar("alA", "alpha*A");    
  params.addRequiredCoupledVar("alrhoA", "alpha*rho*A");
  // Coupled aux variables
  params.addCoupledVar("area", 1., "area");
  // Boolean
  params.addParam<bool>("isLiquid", true,"is the fluid liquid or not");

  return params;
}

DensityAux::DensityAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variables
    _alA(coupledValue("alA")),
    _alrhoA(coupledValue("alrhoA")),
    // Coupled aux variables
    _area(coupledValue("area")),
    // Boolean
    _isLiquid(getParam<bool>("isLiquid"))
{}

Real
DensityAux::computeValue()
{
  // Compute the quantity alpha*A
  Real alA = _isLiquid ? _alA[_qp] : _area[_qp]-_alA[_qp];

  // Return the value of the density:
  return _alrhoA[_qp] / alA;
}
