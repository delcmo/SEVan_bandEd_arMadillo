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

#ifndef ENTROPYVOLUMEFRACTIONAUX_H
#define ENTROPYVOLUMEFRACTIONAUX_H

#include "AuxKernel.h"

class EntropyVolumeFractionAux;

template<>
InputParameters validParams<EntropyVolumeFractionAux>();

class EntropyVolumeFractionAux : public AuxKernel
{
public:

  EntropyVolumeFractionAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  // Coupled aux variable
  VariableValue & _vf_liq;
};

#endif // ENTROPYVOLUMEFRACTIONAUX_H
