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

#include "SbaEnergy.h"
/**
This function computes the convection and source terms of the energy equation for phase k. The phase k is assumed to exchange energy with a phase denoted by the index j.
 **/
template<>
InputParameters validParams<SbaEnergy>()
{
  InputParameters params = validParams<Kernel>();

  // Conservative variables for phase k:
  params.addRequiredCoupledVar("alrhoA_k", " phase k alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x_k", "x component of phase k alpha*rho*u*A");
  params.addRequiredCoupledVar("alrhoEA_k", "phase k alpha*rho*u*A");
  // Conservative variables for phase j:
  params.addRequiredCoupledVar("alrhoA_j", " phase k alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x_j", "x component of phase j alpha*rho*u*A");
  params.addRequiredCoupledVar("alrhoEA_j", "phase j alpha*rho*u*A");
  // Coupled aux variables:
  params.addCoupledVar("area", 1., "cross-section");
  params.addRequiredCoupledVar("liquid_volume_fraction", "liquid volume fraction");
  // Parameters:
  params.addParam<bool>("isLiquid", true, "boolean to determine if liquid phase or not");
  params.addParam<bool>("interfacial_variables_on", false, "boolean for pressure relaxation");  
  params.addParam<RealVectorValue>("gravity", (0., 0., 0.), "gravity vector");
  // Equation of states:
  params.addRequiredParam<UserObjectName>("eos_k", "Equation of state for phase k");
  params.addRequiredParam<UserObjectName>("eos_j", "Equation of state for phase j");

  return params;
}

SbaEnergy::SbaEnergy(const std::string & name,
                       InputParameters parameters) :
  Kernel(name, parameters),
    // Conservative variables for phase k:
    _alrhoA_k(coupledValue("alrhoA_k")),
    _alrhouA_x_k(coupledValue("alrhouA_x_k")),
    _alrhoEA_k(coupledValue("alrhoEA_k")),
    // Conservative variables for phase j:
    _alrhoA_j(coupledValue("alrhoA_j")),
    _alrhouA_x_j(coupledValue("alrhouA_x_j")),
    _alrhoEA_j(coupledValue("alrhoEA_j")),
    // Coupled aux variables
    _area(coupledValue("area")),
    _liquid_vf(coupledValue("liquid_volume_fraction")),
    _grad_liquid_vf(coupledGradient("liquid_volume_fraction")),
    // Parameters:
    _isLiquid(getParam<bool>("isLiquid")),
    _gravity(getParam<RealVectorValue>("gravity")),
    // Equation of states:
    _eos_k(getUserObject<EquationOfState>("eos_k")),
    _eos_j(getUserObject<EquationOfState>("eos_j")),
    // Interfacial variables
    _PI(getMaterialProperty<Real>("interfacial_pressure")),
    _PI_bar(getMaterialProperty<Real>("average_interfacial_pressure")),
    _velI(getMaterialProperty<Real>("interfacial_velocity")),
    _velI_bar(getMaterialProperty<Real>("average_interfacial_velocity")),
    // Relaxation coefficients
    _P_rel(getMaterialProperty<Real>("pressure_relaxation")),
    _vel_rel(getMaterialProperty<Real>("velocity_relaxation"))
{
  mooseAssert(_mesh.dimension() != 1, "The function "<<this->name()<<" can only be used with a 1-D mesh");
}

Real SbaEnergy::computeQpResidual()
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
  Real rho_j = _alrhoA_j[_qp] / (alpha_j*_area[_qp]);

  // Compute total energies:
  Real rhoE_k = _alrhoEA_k[_qp] / (alpha_k*_area[_qp]);
  Real rhoE_j = _alrhoEA_j[_qp] / (alpha_j*_area[_qp]);

  // Compute pressures:
  Real pressure_k = _eos_k.pressure(rho_k, rho_k*vel_k, rhoE_k);
  Real pressure_j = _eos_j.pressure(rho_j, rho_j*vel_j, rhoE_j);

  // Compute convective term:
  Real conv_k = alpha_k*_area[_qp]*vel_k*(rhoE_k+pressure_k);

  // Compute volume fraction and area gradient terms:
  Real grad_term = _PI[_qp]*_velI[_qp];
  grad_term *= _area[_qp]*grad_alpha_k;

  // Velocity relaxation term:
  Real vel_rel_term = _area[_qp]*_velI_bar[_qp]*_vel_rel[_qp]*(vel_j - vel_k);

  // Pressure relaxation term:
  Real press_rel_term = _area[_qp]*_PI_bar[_qp]*_P_rel[_qp]*(pressure_j-pressure_k);

  // Relaxation terms:
  Real rel_term = press_rel_term+vel_rel_term;

  // Gravity work:
  Real gravity_term = _alrhoA_k[_qp]*vel_k*_gravity(0);

  // Return
  return -conv_k*_grad_test[_i][_qp](0) - (rel_term+grad_term+gravity_term)*_test[_i][_qp];
}

Real SbaEnergy::computeQpJacobian()
{
    return 0;
}

Real SbaEnergy::computeQpOffDiagJacobian( unsigned int _jvar)
{
    return 0;
}