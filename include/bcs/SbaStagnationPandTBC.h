#ifndef SBASTAGNATIONPANDTBC_H
#define SBASTAGNATIONPANDTBC_H

#include "IntegratedBC.h"
#include "StiffenedGasEquationOfState.h"

// Forward Declarations
class SbaStagnationPandTBC;
class StiffenedGasEquationOfState;

template<>
InputParameters validParams<SbaStagnationPandTBC>();

/**
 * The boundary condition with specified stagnation pressure and temperature
 */

class SbaStagnationPandTBC : public IntegratedBC
{

public:
  SbaStagnationPandTBC(const std::string & name, InputParameters parameters);

  virtual ~SbaStagnationPandTBC(){}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);

  enum EFlowEquationType
  {
    CONTINUITY = 0,
    XMOMENTUM = 1,
    ENERGY = 2
  };

  // which equation (mass/momentum/energy) this BC is acting on
  MooseEnum _eqn_type;

  // Coupled variables
  VariableValue & _alA;  
  VariableValue & _alrhoA;
  VariableValue & _alrhouA_x;

  // Coupled aux variables:
  VariableValue & _area;

  // Specified stagnation variables:
  Real _p0_bc;
  Real _T0_bc;

  // Calculated rho_0, K, etc. on the boundary:
  Real _rho0_bc;
  Real _H0_bc;
  Real _K;
  Real _H_bar;

  // Equation of state:
  const StiffenedGasEquationOfState & _eos;

  // Boolean phase
  bool _isLiquid;
};

#endif // SbaStagnationPandTBC_H