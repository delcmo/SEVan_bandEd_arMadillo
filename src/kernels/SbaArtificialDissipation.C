/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "SbaArtificialDissipation.h"
/**
This function computes the dissipative terms for all of the 1-d equations.
 */
template<>
InputParameters validParams<SbaArtificialDissipation>()
{
  InputParameters params = validParams<Kernel>();

  // Equation and diffusion names:
  params.addParam<std::string>("equation_name", "INVALID", "Name of the equation.");
  // Coupled aux variables
  params.addRequiredCoupledVar("density", "density of the fluid");
  params.addRequiredCoupledVar("pressure", "pressure of the fluid");
  params.addRequiredCoupledVar("velocity_x", "x component of the velocity");
  params.addRequiredCoupledVar("internal_energy", "internal energy of the fluid");
  params.addCoupledVar("area", 1., "area of the geometry");
  params.addRequiredCoupledVar("liquid_volume_fraction", "liquid volume fraction");
  // Boolean
  params.addParam<bool>("isLiquid", true, "boolean to determine if liquid phase or not");

  return params;
}

SbaArtificialDissipation::SbaArtificialDissipation(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Declare equation types
    _equ_type("VOLUME_FRACTION CONTINUITY XMOMENTUM ENERGY INVALID", getParam<std::string>("equation_name")),
    // Boolean
    _isLiquid(getParam<bool>("isLiquid")),
    // Coupled auxilary variables
    _rho(coupledValue("density")),
    _pressure(coupledValue("pressure")),
    _grad_rho(coupledGradient("density")),
    _grad_press(coupledGradient("pressure")),
    _vel_x(coupledValue("velocity_x")),
    _grad_vel_x(coupledGradient("velocity_x")),
    _rhoe(coupledValue("internal_energy")),
    _grad_rhoe(coupledGradient("internal_energy")),
    _area(coupledValue("area")),
    _alpha_liq(coupledValue("liquid_volume_fraction")),
    _grad_alpha_liq(coupledGradient("liquid_volume_fraction")),
    // Get material property: viscosity coefficient.
    _mu(_isLiquid ? getMaterialProperty<Real>("mu_liq") : getMaterialProperty<Real>("mu_gas")),
    _kappa(_isLiquid ? getMaterialProperty<Real>("kappa_liq") : getMaterialProperty<Real>("kappa_gas")),
    _beta(_isLiquid ? getMaterialProperty<Real>("beta_liq") : getMaterialProperty<Real>("beta_gas"))
{
  mooseAssert(_mesh.dimension() != 1, "The function "<<this->name()<<" can only be used with a 1-D mesh");
}

Real SbaArtificialDissipation::computeQpResidual()
{
  // Initialize the artificial dissipative flux
  Real flux(0.);

  // Determine if cell is on boundary or not:
//  if (_current_elem->node(_i) != 0 || _current_elem->node(_i) != _mesh.nNodes()-1) // THIS IS NOT A ROBUST IMPLEMENTATION
//  {
    // Phase void fraction:
    Real alpha = _isLiquid ? _alpha_liq[_qp] : (1-_alpha_liq[_qp]);
    Real grad_alpha = _isLiquid ? _grad_alpha_liq[_qp](0) : -_grad_alpha_liq[_qp](0);

    // Compute l = beta * grad(alpha):
    Real l_k = grad_alpha;
    l_k *= _beta[_qp];

    // Compute f = kappa * grad(rho):
    Real f_k = _grad_rho[_qp](0);
    f_k *= alpha * _kappa[_qp];
    f_k += _rho[_qp] * l_k;

    // Compute g = mu * grad(vel):
    Real g_k = _grad_vel_x[_qp](0);
    g_k *= _mu[_qp] * _rho[_qp] * alpha;

    // Compute h = kappa * grad(rho*e):
    Real h_k = _grad_rhoe[_qp](0);
    h_k *= alpha * _kappa[_qp];

    // Compute the artificial dissipative flux:
    switch (_equ_type)
    {
      case VOLUME_FRACTION:
        flux = l_k;
        break;
      case CONTINUITY:
        flux = f_k;
        break;
      case XMOMENTUM:
        flux = _vel_x[_qp]*f_k + _mu[_qp]*alpha*_rho[_qp]*_grad_vel_x[_qp](0);
        break;
      case ENERGY:
        flux = h_k + 0.5*f_k*_vel_x[_qp]*_vel_x[_qp] + _vel_x[_qp]*g_k +  _rhoe[_qp]*l_k;
        break;
      default:
        mooseError("INVALID equation name.");
    }
//  }
//  else
//    flux = 0.;

  // Return value
  return _area[_qp]*flux*_grad_test[_i][_qp](0);
}

Real SbaArtificialDissipation::computeQpJacobian()
{
  return 0.;
}

Real SbaArtificialDissipation::computeQpOffDiagJacobian( unsigned int _jvar)
{
  return 0.*_jvar;
}
