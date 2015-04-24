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

#include "SbaMomentum.h"
/**
This function computes the convection and source terms of the momentum equation for phase k. The phase k is assumed to exchange momentum with a phase denoted by the index j.
 **/
template<>
InputParameters validParams<SbaMomentum>()
{
  InputParameters params = validParams<Kernel>();

  // Conservative variables for phase k:
  params.addRequiredCoupledVar("alrhoA_k", " phase k alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x_k", "x component of phase k alpha*rho*u*A");
  params.addRequiredCoupledVar("alrhoEA_k", "phase k alpha*rho*u*A");
  // Conservative variables for phase j:
  params.addRequiredCoupledVar("alrhoA_j", " phase k alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x_j", "x component of phase j alpha*rho*u*A");
  // Coupled aux variables:
  params.addCoupledVar("area", 1., "cross-section");
  params.addRequiredCoupledVar("liquid_volume_fraction", "liquid volume fraction");
  // Parameters:
  params.addParam<bool>("isLiquid", true, "boolean to determine if liquid phase or not");
  params.addParam<RealVectorValue>("gravity", (0., 0., 0.), "gravity vector");
  // Equation of states:
  params.addRequiredParam<UserObjectName>("eos_k", "Equation of state for phase k");

  return params;
}

SbaMomentum::SbaMomentum(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Conservative variables for phase k:
    _alrhoA_k(coupledValue("alrhoA_k")),
    _alrhouA_x_k(coupledValue("alrhouA_x_k")),
    _alrhoEA_k(coupledValue("alrhoEA_k")),
    // Conservative variables for phase j:
    _alrhoA_j(coupledValue("alrhoA_j")),
    _alrhouA_x_j(coupledValue("alrhouA_x_j")),
    // Coupled aux variables
    _area(coupledValue("area")),
    _grad_area(coupledGradient("area")),
    _liquid_vf(coupledValue("liquid_volume_fraction")),
    _grad_liquid_vf(coupledGradient("liquid_volume_fraction")),
    // Parameters:
    _isLiquid(getParam<bool>("isLiquid")),
    _gravity(getParam<RealVectorValue>("gravity")),
    // Equation of states:
    _eos_k(getUserObject<EquationOfState>("eos_k")),
    // Interfacial variables
    _PI(getMaterialProperty<Real>("interfacial_pressure")),
    // Relaxation coefficients
    _vel_rel(getMaterialProperty<Real>("velocity_relaxation"))
{
  if (_mesh.dimension() != 1)
    mooseError("The function "<<this->name()<<" can only be used with a 1-D mesh");
}

Real SbaMomentum::computeQpResidual()
{
  // Compute volume fraction and its derivative of the phase k (liquid or gas):
  Real alpha_k = _isLiquid ? _liquid_vf[_qp] : 1.-_liquid_vf[_qp];
  Real alpha_j = 1.-alpha_k;
  Real grad_alpha_k =_isLiquid ? _grad_liquid_vf[_qp](0) : -_grad_liquid_vf[_qp](0);

  // Compute velocity vectors:
  Real vel_k = _alrhouA_x_k[_qp] / _alrhoA_k[_qp];
  Real vel_j = _alrhouA_x_j[_qp] / _alrhoA_j[_qp];

  // Compute densities:
  Real rho_k = _alrhoA_k[_qp] / (alpha_k*_area[_qp]);

  // Compute pressure for phase k:
  Real rhoE_k = _alrhoEA_k[_qp] / (alpha_k*_area[_qp]);  
  Real pressure_k = _eos_k.pressure(rho_k, vel_k, rhoE_k);

  // Compute convective term:
  Real conv_k = _alrhouA_x_k[_qp]*vel_k+alpha_k*vel_k*pressure_k;

  // Compute volume fraction and area gradient terms:
  Real grad_term = _area[_qp]*_PI[_qp]*grad_alpha_k;
  grad_term += pressure_k*alpha_k*_grad_area[_qp](0);

  // Velocity relaxation term:
  Real vel_rel_term = _area[_qp]*_vel_rel[_qp]*(vel_j - vel_k);

  // Gravity term:
  Real gravity_term = _alrhoA_k[_qp]*_gravity(0);

  // Return
  return -conv_k*_grad_test[_i][_qp](0) - (vel_rel_term+grad_term+gravity_term)*_test[_i][_qp];
}

Real SbaMomentum::computeQpJacobian()
{
    return 0;
}

Real SbaMomentum::computeQpOffDiagJacobian( unsigned int _jvar)
{
    return 0;
}
