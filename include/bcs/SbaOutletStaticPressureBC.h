#ifndef SBAOUTLETSTATICPRESSUREBC_H
#define SBAOUTLETSTATICPRESSUREBC_H

#include "IntegratedBC.h"
#include "StiffenedGasEquationOfState.h"

// Forward Declarations
class SbaOutletStaticPressureBC;
class StiffenedGasEquationOfState;

template<>
InputParameters validParams<SbaOutletStaticPressureBC>();

class SbaOutletStaticPressureBC : public IntegratedBC
{

public:
  SbaOutletStaticPressureBC(const std::string & name, InputParameters parameters);

  virtual ~SbaOutletStaticPressureBC(){}

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
  VariableValue & _alrhoEA;

  // Coupled aux variables
  VariableValue & _area;

  // Specified pressure and temperature values
  Real _p_bc;

  // Equation of state
  const StiffenedGasEquationOfState & _eos;

  // Boolean phase
  bool _isLiquid;
};

#endif // SBAOUTLETSTATICPRESSUREBC_H

