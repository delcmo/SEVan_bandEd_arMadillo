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

#ifndef SBAARTIFICIALDISSIPATION_H
#define SBAARTIFICIALDISSIPATION_H

#include "Kernel.h"

// Forward Declarations
class SbaArtificialDissipation;

template<>
InputParameters validParams<SbaArtificialDissipation>();

class SbaArtificialDissipation : public Kernel
{
public:

  SbaArtificialDissipation(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int _jvar);
    
private:
  // Equations types
  enum EquationType
  {
      VOLUME_FRACTION = 0,
      CONTINUITY = 1,
      XMOMENTUM = 2,
      ENERGY = 3
  };

  // Diffusion name
  MooseEnum _equ_type;

  // Boolean for phase
  bool _isLiquid;

  // Coupled aux variables:
  VariableValue & _rho;
  VariableValue & _pressure;
  VariableGradient & _grad_rho;
  VariableGradient & _grad_press;
  VariableValue & _vel_x;
  VariableGradient & _grad_vel_x;
  VariableValue & _rhoe;
  VariableGradient & _grad_rhoe;
  VariableValue & _area;
  VariableValue & _alpha_liq;
  VariableGradient & _grad_alpha_liq;

  // Material property: viscosity coefficient.
  const MaterialProperty<Real> & _mu;
  const MaterialProperty<Real> & _kappa;
  const MaterialProperty<Real> & _beta;
};

#endif // SBAARTIFICIALDISSIPATION_H
