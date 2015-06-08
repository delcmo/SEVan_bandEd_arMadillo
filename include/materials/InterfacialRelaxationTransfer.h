#ifndef INTERFACIALRELAXATIONTRANSFER_H
#define INTERFACIALRELAXATIONTRANSFER_H

#include "Material.h"
#include "MaterialProperty.h"
#include "EquationOfState.h"

//Forward Declarations
class InterfacialRelaxationTransfer;

template<>
InputParameters validParams<InterfacialRelaxationTransfer>();

class InterfacialRelaxationTransfer : public Material
{
public:
  InterfacialRelaxationTransfer(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpProperties();

private:
    // Variable for definition of the interfacial variables:
    enum InterfacialType
    {
      BERRY = 0,
      AMBROSSO = 1,
      LIANG = 2,
      NO_RELAXATION = 3
    };
    MooseEnum _interfacial_param_def;
    
    Real _xi;
    
  // Coupled variables phase k
  VariableValue & _alrhoA_k;
  VariableValue & _alrhouA_x_k;
  VariableValue & _alrhoEA_k;

  // Coupled variables phase j
  VariableValue & _alrhoA_j;
  VariableValue & _alrhouA_x_j;
  VariableValue & _alrhoEA_j;

  // Coupled aux variables:
  VariableValue & _vf_k;
  VariableGradient & _grad_vf_k;
  VariableValue & _area;

  // Interfacial variables:
  MaterialProperty<Real> & _Aint;
  MaterialProperty<Real> & _PI;
  MaterialProperty<Real> & _velI;
  MaterialProperty<Real> & _PI_bar;
  MaterialProperty<Real> & _velI_bar;

  // Relaxation parameters:
  MaterialProperty<Real> & _press_rel_coeff;
  MaterialProperty<Real> & _vel_rel_coeff;

  // Maximum specific interfacial area:
  Real _Aint_max_press;
  Real _Aint_max_vel;

  // Equation of state:
  const EquationOfState & _eos_k;
  const EquationOfState & _eos_j;

  // Booleans
  bool _has_max_specific_interf_area;
};

#endif //INTERFACIALRELAXATIONTRANSFER_H