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

#ifndef SBAMOMENTUM_H
#define SBAMOMENTUM_H

#include "Kernel.h"
#include "EquationOfState.h"

// Forward Declarations
class SbaMomentum;

template<>
InputParameters validParams<SbaMomentum>();

class SbaMomentum : public Kernel
{
public:

  SbaMomentum(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int _jvar);

private:
  // Coupled variables phase k
  VariableValue & _alrhoA_k;
  VariableValue & _alrhouA_x_k;
  VariableValue & _alrhoEA_k;

  // Coupled variables phase j
  VariableValue & _alrhoA_j;
  VariableValue & _alrhouA_x_j;
    
  // Coupled aux variables:
  VariableValue & _area;
  VariableGradient & _grad_area;
  VariableValue & _liquid_vf;
  VariableGradient & _grad_liquid_vf;

  // Parameters
  bool _isLiquid;
  RealVectorValue _gravity;

  // Equation of state:
  const EquationOfState & _eos_k;
    
  // Interfacial variables
  MaterialProperty<Real> & _PI;

  // Relaxation parameters.
  MaterialProperty<Real> & _vel_rel;
};

#endif // SBAMOMENTUM_H