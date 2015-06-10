#ifndef EQUATIONOFSTATE_H
#define EQUATIONOFSTATE_H

#include "GeneralUserObject.h"

// Forward Declarations
class EquationOfState;

template<>
InputParameters validParams<EquationOfState>();

class EquationOfState : public GeneralUserObject
{
public:
  // Constructor
  EquationOfState(const std::string & name, InputParameters parameters);

  // Destructor
  virtual ~EquationOfState();

  /**
   * Called when this object needs to compute something.
   */
  virtual void execute() {};

  /**
   * Called before execute() is ever called so that data can be cleared.
   */
  virtual void initialize(){};

  /**
   * Finalize.  This is called _after_ execute() and _after_ threadJoin()!  This is probably where you want to do MPI communication!
   */
  virtual void finalize() {};

  // The interface for derived EquationOfState objects to implement...
  virtual Real pressure(Real rho, Real vel, Real rhoE) const;

  // The interface for derived EquationOfState objects to implement...
  virtual Real temperature(Real rho, Real vel, Real rhoE) const;

  // Sound speed squared
  virtual Real c2(Real rho, Real vel, Real rhoE) const;

  // Sound speed squared
  virtual Real c2_from_rho_p(Real rho, Real pressure) const;

  // Computes density from pressure and temperature for single-phase
  virtual Real rho_from_p_T(Real pressure, Real temperature) const;

  // Computes pressure from density and temperature for single-phase  
  virtual Real p_from_rho_T(Real rho, Real temperature) const;

  // Compute internal energy from pressure and density
  virtual Real e_from_p_rho(Real pressure, Real rho) const;

  // The interface for derived EquationOfState objects to implement...
  virtual Real temp_from_p_rho(Real pressure, Real rho) const;

protected:
  // Prints an error message for non-implemented functions
  void error_not_implemented(std::string method_name) const;
};

#endif // EQUATIONOFSTATE_H