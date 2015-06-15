#ifndef SBADIRICHLETBC_H
#define SBADIRICHLETBC_H

#include "NodalBC.h"
#include "EquationOfState.h"

//Forward Declarations
class SbaDirichletBC;

template<>
InputParameters validParams<SbaDirichletBC>();

/**
 * Implements space-dependent Dirichlet BC.
 */
class SbaDirichletBC : public NodalBC
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  SbaDirichletBC(const std::string & name, InputParameters parameters);

  virtual ~SbaDirichletBC(){}

protected:

  virtual Real computeQpResidual();

    enum EFlowEquationType
    {
        VOIDFRACTION = 0,
        CONTINUITY = 1,
        XMOMENTUM = 2,
        ENERGY = 3
    };

  // Initial condition type
  int _ics_type;

  // Function area
  Function * _area;

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

  // Equation of state:
  const EquationOfState & _eos;

  // Booleans:
  bool _isLiquid;
  bool _isArea;
  bool _isLeft;

};

#endif // SBADIRICHLETBC_H