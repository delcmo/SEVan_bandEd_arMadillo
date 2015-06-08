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
This function computes the void fraction from the variables 'alhpaA' and 'A'.
**/

#include "VolumeFractionAux.h"

template<>
InputParameters validParams<VolumeFractionAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Coupled variable
  params.addRequiredCoupledVar("alA", "alpha*A");
  // Coupled aux variable
  params.addCoupledVar("area", 1., "area");

  return params;
}

VolumeFractionAux::VolumeFractionAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Coupled variable
    _alA(coupledValue("alA")),
    // Coupled aux variable
    _area(coupledValue("area"))
{}

Real
VolumeFractionAux::computeValue()
{
  // Compute the volume fraction
  Real alpha = _alA[_qp] / _area[_qp];

  // return the value
  if ( alpha < 0 )
      return 0.;
  else if ( alpha > 1 )
      return 1.;
  else
      return alpha;
}
