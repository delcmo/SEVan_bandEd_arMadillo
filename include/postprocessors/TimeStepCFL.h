#ifndef TIMESTEPCFL_H
#define TIMESTEPCFL_H

#include "ElementPostprocessor.h"
#include "EquationOfState.h"

class TimeStepCFL;

template<>
InputParameters validParams<TimeStepCFL>();

/**
 * The inviscid time step stability limit:
 *
 * h_e \over {|\vec u| + c}
 */
class TimeStepCFL : public ElementPostprocessor
{
public:
  TimeStepCFL(const std::string & name, InputParameters parameters);
  virtual ~TimeStepCFL();

  virtual void initialize();
  virtual void execute();
  virtual Real getValue();
  virtual void threadJoin(const UserObject & uo);

protected:

  // Coupled variables
  VariableValue & _alA_k;
  VariableValue & _alrhoA_k;
  VariableValue & _alrhouA_x_k;
  VariableValue & _alrhoEA_k;
  VariableValue & _alrhoA_j;
  VariableValue & _alrhouA_x_j;
  VariableValue & _alrhoEA_j;

  // Coupled aux variables
  VariableValue & _area;

  // Equation of state
  const EquationOfState & _eos_k;
  const EquationOfState & _eos_j;

  // Material property
  const MaterialProperty<Real> & _velI;

  // Parameter
  Real _cfl;
  Real _value;
};


#endif // TIMESTEPCFL_H
