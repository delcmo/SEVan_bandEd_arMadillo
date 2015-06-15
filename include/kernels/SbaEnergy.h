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

#ifndef SBAENERGY_H
#define SBAENERGY_H

#include "Kernel.h"
#include "EquationOfState.h"

// Forward Declarations
class SbaEnergy;

template<>
InputParameters validParams<SbaEnergy>();

class SbaEnergy : public Kernel
{
public:

  SbaEnergy(const std::string & name,
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
  VariableValue & _alrhoEA_j;
    
  // Coupled aux variables:
  VariableValue & _area;
  VariableValue & _liquid_vf;
  VariableGradient & _grad_liquid_vf;

  // Parameters
  bool _isLiquid;
  RealVectorValue _gravity;

  // Equation of state:
  const EquationOfState & _eos_k;
  const EquationOfState & _eos_j;
    
  // Interfacial variables
  const MaterialProperty<Real> & _PI;
  const MaterialProperty<Real> & _PI_bar;
  const MaterialProperty<Real> & _velI;
  const MaterialProperty<Real> & _velI_bar;

  // Relaxation parameters.
  const MaterialProperty<Real> & _P_rel;
  const MaterialProperty<Real> & _vel_rel;
};

#endif // SBAENERGY_H