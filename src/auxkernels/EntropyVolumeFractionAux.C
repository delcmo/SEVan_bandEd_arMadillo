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
Computes the entropy volume fraction.
**/
#include "EntropyVolumeFractionAux.h"

template<>
InputParameters validParams<EntropyVolumeFractionAux>()
{
  InputParameters params = validParams<AuxKernel>();

  // Aux variables:
  params.addRequiredCoupledVar("volume_fraction_liquid", "liquid volume fraction");

  return params;
}

EntropyVolumeFractionAux::EntropyVolumeFractionAux(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    // Aux variables:
    _vf_liq(coupledValue("volume_fraction_liquid"))

{
}

Real
EntropyVolumeFractionAux::computeValue()
{
  // Return the entropy value
  return 0.5*_vf_liq[_qp]*_vf_liq[_qp];
}
