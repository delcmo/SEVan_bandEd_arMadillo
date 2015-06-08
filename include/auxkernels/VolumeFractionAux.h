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

#ifndef VOLUMEFRACTIONAUX_H
#define VOLUMEFRACTIONAUX_H

#include "AuxKernel.h"


//Forward Declarations
class VolumeFractionAux;

template<>
InputParameters validParams<VolumeFractionAux>();

/**
 * Coupled auxiliary value
 */
class VolumeFractionAux : public AuxKernel
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  VolumeFractionAux(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();

  // Coupled variable
  VariableValue & _alA;

  // Coupled aux variable
  VariableValue & _area;
};

#endif // VOLUMEFRACTIONAUX_H
