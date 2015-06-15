#include "SbaStagnationPandTBC.h"

template<>
InputParameters validParams<SbaStagnationPandTBC>()
{
  InputParameters params = validParams<IntegratedBC>();

  params.addRequiredParam<std::string>("equation_name", "The name of the equation this BC is acting on");
  // Coupled conservative variables:
  params.addCoupledVar("alA", "alpha*A");  
  params.addCoupledVar("alrhoA", "alpha*rho*A");
  params.addCoupledVar("alrhouA_x", "");
  // Coupled aux variables:
  params.addCoupledVar("area", "Coupled area variable");
  // Pressure and temperature stagnation values
  params.addRequiredParam<Real>("p0_bc", "Stagnation pressure at the boundary");
  params.addRequiredParam<Real>("T0_bc", "Liquid stagnation temperature at the boundary");
  // Equation of state
  params.addRequiredParam<UserObjectName>("eos", "The name of equation of state object to use.");
  // Boolean
  params.addParam<bool>("isLiquid", true, "is liquid or not?");

  return params;
}

SbaStagnationPandTBC::SbaStagnationPandTBC(const std::string & name, InputParameters parameters) :
    IntegratedBC(name, parameters),
    // Type of equation:
    _eqn_type("CONTINUITY XMOMENTUM ENERGY INVALID", getParam<std::string>("equation_name")),
    // Coupled aux variables:
    _alA(coupledValue("alA")),
    _alrhoA(coupledValue("alrhoA")),
    _alrhouA_x(coupledValue("alrhouA_x")),
    // Coupled aux variables:
    _area(coupledValue("area")),
    // Stagnation variables:
    _p0_bc(getParam<Real>("p0_bc")),
    _T0_bc(getParam<Real>("T0_bc")),
    // Equation of state:
    _eos(getUserObject<StiffenedGasEquationOfState>("eos")),
    // Boolean:
    _isLiquid(getParam<bool>("isLiquid"))
{
  mooseAssert(_mesh.dimension() != 1, "The function "<<this->name()<<" can only be used with a 1-D mesh");

  // Pre-compute some stagnation cpefficients:
  _rho0_bc = _eos.rho_from_p_T(_p0_bc, _T0_bc);
  _H0_bc = _eos.e_from_p_rho(_p0_bc, _rho0_bc) + _p0_bc / _rho0_bc;
  _K = (_p0_bc + _eos.pInf()) / std::pow(_rho0_bc, _eos.gamma());
  _H_bar = _eos.gamma() * (_p0_bc + _eos.pInf()) / _rho0_bc / (_eos.gamma() - 1);
}

Real
SbaStagnationPandTBC::computeQpResidual()
{
  // Compute the void fraction:
  Real alpha_k = _isLiquid ? _alA[_qp]/_area[_qp] : 1-_alA[_qp]/_area[_qp];

  // Compute u_star and v_star:
  Real vel_k = _alrhouA_x[_qp]/_alrhoA[_qp];

  // Compute rho_star and static pressure:
  Real rho_k_star = std::pow((_H_bar - 0.5*vel_k*vel_k)*(_eos.gamma()-1)/(_eos.gamma())/_K, 1./(_eos.gamma()-1));
  Real p_k_bc = _K * std::pow(rho_k_star, _eos.gamma()) - _eos.pInf();

  // return value
  switch (_eqn_type)
  {
    case CONTINUITY:
      return alpha_k*rho_k_star*vel_k*_area[_qp]*_normals[_qp](0)*_test[_i][_qp];
      break;
    case XMOMENTUM:
      return alpha_k*_area[_qp]*(rho_k_star*vel_k*vel_k+p_k_bc)*_normals[_qp](0)*_test[_i][_qp];
      break;
    case ENERGY:
      return alpha_k*rho_k_star*vel_k*_area[_qp]*_H0_bc*_normals[_qp](0)*_test[_i][_qp];
      break;      
    default:
      mooseError("The equation name given in the input file is not supported in the \"SbaStagnationPandTBC\" type of boundary condition.");
      return 0.;
  }
}

Real
SbaStagnationPandTBC::computeQpJacobian()
{
  // TODO
  return 0;
}

Real
SbaStagnationPandTBC::computeQpOffDiagJacobian(unsigned jvar)
{
  // TODO
  return 0;
}
