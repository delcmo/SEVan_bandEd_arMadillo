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

#include "SbaVolumeFraction.h"
/**
This function computes the convection and source terms of the volume fraction equation for phase k. The phase k is assumed to exchange volume fraction with a phase denoted by the index j.
 **/
template<>
InputParameters validParams<SbaVolumeFraction>()
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
  // Equation of states:
  params.addRequiredParam<UserObjectName>("eos_k", "Equation of state for phase k");
  params.addRequiredParam<UserObjectName>("eos_j", "Equation of state for phase j");

  return params;
}

SbaVolumeFraction::SbaVolumeFraction(const std::string & name,
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
    // Equation of states:
    _eos_k(getUserObject<EquationOfState>("eos_k")),
    _eos_j(getUserObject<EquationOfState>("eos_j")),
    // Interfacial variables
    _velI(getMaterialProperty<Real>("interfacial_velocity")),
    // Relaxation coefficients
    _P_rel(getMaterialProperty<Real>("pressure_relaxation"))
{
  if (_mesh.dimension() != 1)
    mooseError("The function"<<this->name()<<" can only be used with a 1-D mesh");
}

Real SbaVolumeFraction::computeQpResidual()
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
  Real pressure_k = _eos_k.pressure(rho_k, vel_k, rhoE_k);
  Real pressure_j = _eos_j.pressure(rho_j, vel_j, rhoE_j);

  // Compute convective term:
  Real conv_k = _area[_qp]*_velI[_qp]*grad_alpha_k;

  // Pressure relaxation term:
  Real press_rel_term = _area[_qp]*_P_rel[_qp]*(pressure_j-pressure_k);

  // Return
  return -conv_k*_test[_i][_qp] - press_rel_term*_test[_i][_qp];
}

Real SbaVolumeFraction::computeQpJacobian()
{
    return 0;
}

Real SbaVolumeFraction::computeQpOffDiagJacobian( unsigned int _jvar)
{
    return 0;
}
