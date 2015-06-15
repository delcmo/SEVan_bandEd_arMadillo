#ifndef COMPUTEVISCOSITYCOEFFICIENT_H
#define COMPUTEVISCOSITYCOEFFICIENT_H

#include "Material.h"
#include "MaterialProperty.h"
#include "EquationOfState.h"

//Forward Declarations
class ComputeViscosityCoefficient;

template<>
InputParameters validParams<ComputeViscosityCoefficient>();

class ComputeViscosityCoefficient : public Material
{
public:
  ComputeViscosityCoefficient(const std::string & name, InputParameters parameters);
  virtual ~ComputeViscosityCoefficient();

protected:
  virtual void initQpStatefulProperties();  
  virtual void computeQpProperties();

private:    
  // Boolean
  bool _is_liquid;
  bool _is_jump_on;
  bool _is_first_order_visc;

  // Conservative variables
  VariableValue & _alrhoA_k;
  VariableValue & _alrhouA_x_k;

  // Liquid void fraction:
  VariableValue & _ent_vf_liq;
  VariableValue & _ent_vf_liq_old;
  VariableValue & _ent_vf_liq_older;
  VariableGradient & _grad_ent_vf_liq;

  // Pressure:
  VariableValue & _press;
  VariableValue & _press_old;
  VariableValue & _press_older;
  VariableGradient & _grad_press;

  // Density:
  VariableValue & _rho;
  VariableValue & _rho_old;
  VariableValue & _rho_older;
  VariableGradient & _grad_rho;

  // Variables for jump:
  VariableValue & _jump_grad_press;
  VariableValue & _jump_grad_dens;
  VariableValue & _jump_grad_vf;

  // Material property: interfacial velocity.
  MaterialProperty<Real> & _velI;

  // Multiplicative coefficient for viscosity:
  Real _Cmax;
  Real _Ce;
  Real _Cjump;
  Real _Ce_vf;
  Real _Cjump_vf;

  // UserObject: equation of state
  const EquationOfState & _eos;

  // Postprocessors
  std::string _vf_pps_name;

  // Declare material properties: viscosity coefficients.
  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> & _mu_max;
  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _kappa_max;
  MaterialProperty<Real> & _beta;
  MaterialProperty<Real> & _beta_max;
};

#endif // COMPUTEVISCOSITYCOEFFICIENT_H
