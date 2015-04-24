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

#ifndef SBAMASS_H
#define SBAMASS_H

#include "Kernel.h"

// Forward Declarations
class SbaMass;

template<>
InputParameters validParams<SbaMass>();

class SbaMass : public Kernel
{
public:

  SbaMass(const std::string & name,
             InputParameters parameters);

protected:

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian( unsigned int _jvar);

private:
  // Coupled variables phase k
  VariableValue & _alrhouA_x_k;
};

#endif // SBAMASS_H