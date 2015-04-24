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

#ifndef SBAICS_H
#define SBAICS_H

// MOOSE Includes
#include "InitialCondition.h"
#include "EquationOfState.h"
#include "Function.h"

// Forward Declarations
class SbaICs;

template<>
InputParameters validParams<SbaICs>();

/**
 * SbaICs just returns a constant value.
 */
class SbaICs : public InitialCondition
{
public:

  SbaICs(const std::string & name,
            InputParameters parameters);

  virtual Real value(const Point & p);

private:

  // Initial condition type
  int _ics_type;

  // Function area
  Function & _area;

  // Initial values for pressure:
  Real _p_left;
  Real _p_right;

  // Initial values for velocity:
  Real _v_left;
  Real _v_right;

  // Initial values for temperature:
  Real _t_left;
  Real _t_right;

  // Initial values for density:
  Real _rho_left;
  Real _rho_right;

  // Initial values for volume fraction:
  Real _liq_vf_left;
  Real _liq_vf_right;

  // Position of the membrane:
  Real _membrane;
  Real _length;

  // Equation of state:
  const EquationOfState & _eos;

  // Booleans:
  bool _isLiquid;
  bool _isArea;
};

#endif // SBAICS_H