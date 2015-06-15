#include "SbaOutletStaticPressureBC.h"

template<>
InputParameters validParams<SbaOutletStaticPressureBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");

  // Coupled variables:
  params.addRequiredCoupledVar("alA", "alpha*A");
  params.addRequiredCoupledVar("alrhoA", "alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x", "x component of alpha*rho*u*A");
  params.addRequiredCoupledVar("alrhoEA", "alpha*rho*E*A");
  // Coupled aux variables
  params.addCoupledVar("area", 1., "area aux variable");
  // Input parameters:
  params.addRequiredParam<Real>("p_bc", "Static pressure at the boundary");
  // Equation of state:
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  // Boolean
  params.addParam<bool>("isLiquid", true, "is liquid or not?");

  return params;
}

SbaOutletStaticPressureBC::SbaOutletStaticPressureBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Name of the equation:
    _eqn_type("CONTINUITY XMOMENTUM ENERGY INVALID", getParam<std::string>("equation_name")),
    // Coupled variables:
    _alA(coupledValue("alA")),
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA_x(coupledValue("alrhouA_x")),
    _alrhoEA(coupledValue("alrhoEA")),
    // Coupled aux variables
    _area(coupledValue("area")),
    // Boundary condition parameters:
    _p_bc(getParam<Real>("p_bc")),
    // Equation of state:
    _eos(getUserObject<StiffenedGasEquationOfState>("eos")),
    // Boolean:
    _isLiquid(getParam<bool>("isLiquid"))
{

}

Real
SbaOutletStaticPressureBC::computeQpResidual()
{
  // Compute the velocity at the boundary
  Real vel_k_bc = _alrhouA_x[_qp]/_alrhoA[_qp];
  mooseAssert(vel_k_bc*_normals[_qp](0)<0, "The boundary is no longer an outlet.");

  // Compute the void fraction:
  Real alpha_k_bc = _isLiquid ? _alA[_qp]/_area[_qp] : 1.-_alA[_qp]/_area[_qp];

  // Return value
  switch (_eqn_type)
  {
  case CONTINUITY:
    return _alrhouA_x[_qp]*_normals[_qp](0)*_test[_i][_qp];
    break;
  case XMOMENTUM:
    return (_alrhouA_x[_qp]*vel_k_bc+alpha_k_bc*_area[_qp]*_p_bc)*_normals[_qp](0)*_test[_i][_qp];
    break;
  case ENERGY:
    return vel_k_bc*(_alrhoEA[_qp]+alpha_k_bc*_area[_qp]*_p_bc)*_normals[_qp](0)*_test[_i][_qp];
    break;
  default:
    mooseError("The equation name supplied in the input file is not supported in the '"<<name()<<"' type of boundary condition.");
    break;
  }
}

Real
SbaOutletStaticPressureBC::computeQpJacobian()
{
  // TODO
  return 0;
}

Real
SbaOutletStaticPressureBC::computeQpOffDiagJacobian(unsigned jvar)
{
  // TODO
  return 0;
}
