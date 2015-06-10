#ifndef STIFFENEDGASEQUATIONOFSTATE_H
#define STIFFENEDGASEQUATIONOFSTATE_H

#include "EquationOfState.h"

// Forward Declarations
class StiffenedGasEquationOfState;

template<>
InputParameters validParams<StiffenedGasEquationOfState>();


class StiffenedGasEquationOfState : public EquationOfState
{
public:
  StiffenedGasEquationOfState(const std::string & name, InputParameters parameters);
  virtual ~StiffenedGasEquationOfState();

  virtual Real pressure(Real rho, Real rhou, Real rhoE) const;
  virtual Real temperature(Real rho, Real rhou, Real rhoE) const;
  virtual Real c2(Real rho, Real rhou, Real rhoE) const;
  virtual Real c2_from_rho_p(Real, Real) const;
  virtual Real rho_from_p_T(Real pressure, Real temperature) const;
  virtual Real p_from_rho_T(Real rho, Real temperature) const;
  virtual Real e_from_p_rho(Real pressure, Real rho) const;
  virtual Real temp_from_p_rho(Real pressure, Real rho) const;

  virtual Real gamma() const { return _gamma; }
  Real cv() const { return _cv; }
  Real qcoeff() const { return _q; }
  Real pInf() const { return _p_inf; }

protected:
  Real _gamma;
  Real _cv;
  Real _q;
  Real _q_prime;
  Real _p_inf;

  Real _cp;
};

#endif /* STIFFENEDGASEQUATIONOFSTATE_H */
