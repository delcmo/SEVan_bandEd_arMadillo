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

#include "SbaMass.h"
/**
This function computes the convection and source terms of the mass equation for phase k. The phase k is assumed to exchange mass with a phase denoted by the index j.
 **/
template<>
InputParameters validParams<SbaMass>()
{
  InputParameters params = validParams<Kernel>();

  // Conservative variables for phase k:
  params.addRequiredCoupledVar("alrhouA_x_k", "x component of phase k alpha*rho*u*A");

  return params;
}

SbaMass::SbaMass(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Conservative variables for phase k:
    _alrhouA_x_k(coupledValue("alrhouA_x_k"))
{
  if (_mesh.dimension() != 1)
    mooseError("The function "<<this->name()<<" can only be used with a 1-D mesh");
}

Real SbaMass::computeQpResidual()
{
  // Compute convective term:
  Real conv_k = _alrhouA_x_k[_qp];

  // Return
  return -conv_k*_grad_test[_i][_qp](0);
}

Real SbaMass::computeQpJacobian()
{
    return 0;
}

Real SbaMass::computeQpOffDiagJacobian( unsigned int _jvar)
{
    return 0;
}
