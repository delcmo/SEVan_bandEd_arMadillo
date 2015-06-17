#include "ComputeViscosityCoefficient.h"

template<>
InputParameters validParams<ComputeViscosityCoefficient>()
{
  InputParameters params = validParams<Material>();

  // Booleans:
  params.addParam<bool>("is_liquid", true, "the phase is liquid or not.");
  params.addParam<bool>("is_jump_on", true, "use jump or gradients.");
  params.addParam<bool>("is_first_order_visc", true, "use the first-order viscosity coefficients if true.");  
  // Coupled variables
  params.addRequiredCoupledVar("alrhoA_k", " phase k alpha*rho*A");
  params.addRequiredCoupledVar("alrhouA_x_k", "x component of phase k alpha*rho*u*A");
  // Coupled Aux variables:
  params.addRequiredCoupledVar("pressure", "pressure of the fluid");
  params.addRequiredCoupledVar("density", "density of the fluid: rho");
  params.addCoupledVar("entropy_vf_liquid", 1., "liquid void fraction.");
  // Jumps:
  params.addCoupledVar("jump_grad_press", "jump of pressure gradient");
  params.addCoupledVar("jump_grad_rho", "jump of density gradient");
  params.addCoupledVar("jump_grad_vf", "jump of volume fraction gradient");
  // Constant parameter:
  params.addParam<Real>("Cmax", 0.5, "Coefficient for first-order viscosity");
  params.addParam<Real>("Ce", 1., "Coefficient for entropy residual");
  params.addParam<Real>("Cjump", 1., "Coefficient for jump");
  params.addParam<Real>("Ce_vf", 1., "Coefficient for entropy residual of volume fraction equation");
  params.addParam<Real>("Cjump_vf", 1., "Coefficient for jump in volume fraction viscosity coefficient");
  // Userobject:
  params.addRequiredParam<UserObjectName>("eos", "Equation of state");
  // PPS names:
  params.addParam<std::string>("vf_pps_name", "name of the pps for volume fraction viscosity coefficient");

  return params;
}

ComputeViscosityCoefficient::ComputeViscosityCoefficient(const std::string & name, InputParameters parameters) :
    Material(name, parameters),
    // Boolean for phase:
    _is_liquid(getParam<bool>("is_liquid")),
    _is_jump_on(getParam<bool>("is_jump_on")),
    _is_first_order_visc(getParam<bool>("is_first_order_visc")),
    // Conservative variables for phase k:
    _alrhoA_k(coupledValue("alrhoA_k")),
    _alrhouA_x_k(coupledValue("alrhouA_x_k")),
    // Liquid void fraction:
    _ent_vf_liq(_is_liquid ? coupledValue("entropy_vf_liquid") : _zero),
    _ent_vf_liq_old(_is_liquid ? coupledValueOld("entropy_vf_liquid") : _zero),
    _ent_vf_liq_older(_is_liquid ? coupledValueOlder("entropy_vf_liquid") : _zero),
    _grad_ent_vf_liq(_is_liquid ? coupledGradient("entropy_vf_liquid") : _grad_zero),
    // Pressure:
    _press(coupledValue("pressure")),
    _press_old(coupledValueOld("pressure")),
    _press_older(coupledValueOlder("pressure")),
    _grad_press(coupledGradient("pressure")),
    // Density:
    _rho(coupledValue("density")),
    _rho_old(coupledValueOld("density")),
    _rho_older(coupledValueOlder("density")),
    _grad_rho(coupledGradient("density")),
    // Jump of pressure, density and volume fraction gradients:
    _jump_grad_press(isCoupled("jump_grad_press") ? coupledValue("jump_grad_press") : _zero),
    _jump_grad_dens(isCoupled("jump_grad_rho") ? coupledValue("jump_grad_rho") : _zero),
    _jump_grad_vf(isCoupled("jump_grad_vf") ? coupledValue("jump_grad_vf") : _zero),
    // Get interfacial velocity
    _velI(getMaterialProperty<Real>("interfacial_velocity")),
    // Get parameters: Cmax, Ce, Cjump:
    _Cmax(getParam<Real>("Cmax")),
    _Ce(getParam<Real>("Ce")),
    _Cjump(getParam<Real>("Cjump")),
    _Ce_vf(getParam<Real>("Ce_vf")),
    _Cjump_vf(getParam<Real>("Cjump_vf")),
    // UserObject:
    _eos(getUserObject<EquationOfState>("eos")),
    // PPS name:
    _vf_pps_name(getParam<std::string>("vf_pps_name")),
    // Declare material properties used in mass, momentum and energy equations:
    _mu(_is_liquid ? declareProperty<Real>("mu_liq") : declareProperty<Real>("mu_gas")),
    _mu_max(_is_liquid ? declareProperty<Real>("mu_max_liq") : declareProperty<Real>("mu_max_gas")),
    _kappa(_is_liquid ? declareProperty<Real>("kappa_liq") : declareProperty<Real>("kappa_gas")),
    _visc_max(_is_liquid ? declareProperty<Real>("visc_max_liq") : declareProperty<Real>("visc_max_gas")),
    // Declare material property used in volume fraction equation:
    _beta(_is_liquid ? declareProperty<Real>("beta_liq") : declareProperty<Real>("beta_gas")),
    _beta_max(_is_liquid ? declareProperty<Real>("beta_max_liq") : declareProperty<Real>("beta_max_gas"))
{
}

ComputeViscosityCoefficient::~ComputeViscosityCoefficient()
{
}

void
ComputeViscosityCoefficient::initQpStatefulProperties()
{
}

void
ComputeViscosityCoefficient::computeQpProperties()
{
  // Determine h (characteristic length of the current cel):
  Real h = std::pow(_current_elem->volume(), 1./_mesh.dimension());

  // Compute the first-order viscosity coefficients:
  Real c = std::sqrt(_eos.c2_from_rho_p(_rho[_qp], _press[_qp]));
  Real vel = _alrhouA_x_k[_qp]/_alrhoA_k[_qp];
  _mu_max[_qp] = 0.5*h*std::fabs(vel);
  _visc_max[_qp] = 0.5*h*(std::fabs(vel) + c);
  _beta_max[_qp] = 0.5*h*std::fabs(_velI[_qp]);

  // Compute the weights for BDF2 temporal integrator
  Real w0(0.), w1(0.), w2(0.);
  if (_t_step > 1)
  {
    w0 = (2.*_dt+_dt_old)/(_dt*(_dt+_dt_old));
    w1 = -(_dt+_dt_old)/(_dt*_dt_old);
    w2 = _dt/(_dt_old*(_dt+_dt_old));
  }
  else
  {
    w0 =  1. / _dt; w1 = -1. / _dt; w2 = 0.;
  }

  // Compute the Mach number
  Real Mach = std::min(vel/c, 1.);

  // Compute the entropy residual ('res') for the volume fraction equation
  Real res = w0*_ent_vf_liq[_qp]+w1*_ent_vf_liq_old[_qp]+w2*_ent_vf_liq_older[_qp];
  res += _velI[_qp]*_grad_ent_vf_liq[_qp](0);
  res *= _Ce_vf;

  // Compute the jump value
  Real jump = _is_jump_on ? _jump_grad_vf[_qp] : _grad_ent_vf_liq[_qp](0);
  jump *= _Cjump_vf*_velI[_qp];

  // Compute the norm
  Real norm = getPostprocessorValueByName(_vf_pps_name);

  // Compute entropy visc. coeff. for vol. fraction equ.
  Real beta_e = std::max(res, jump);
  beta_e *= h*h/norm;
  _beta_max[_qp] = _is_first_order_visc || _t_step==1 ? _beta_max[_qp] : std::min(beta_e, _beta_max[_qp]);

  // Compute the entropy residual used in the definition of 'kappa' and 'mu': res=DP/Dt-c*c*DrhoDt
  res = -w0*_rho[_qp]-w1*_rho_old[_qp]-w2*_rho_older[_qp];
  res -= vel*_grad_rho[_qp](0);
  res *= c*c;
  res += w0*_press[_qp]+w1*_press_old[_qp]+w2*_press_older[_qp];
  res += vel*_grad_press[_qp](0);

  // Compute the jump values
  Real jump_rho = _is_jump_on ? _jump_grad_dens[_qp] : _grad_rho[_qp](0);
  jump_rho *= c*c;
  Real jump_press = _is_jump_on ? _jump_grad_press[_qp] : _grad_press[_qp](0);
  jump = std::max(jump_rho, jump_press);

  // Compute 'mu'
  norm = Mach*c*c+(1.-Mach)*vel*vel;
  norm *= _rho[_qp];
  Real mu_e = std::max(res, jump);
  mu_e *= h*h/norm;
  _mu[_qp] = _is_first_order_visc || _t_step==1 ? _visc_max[_qp] : std::min(_visc_max[_qp], mu_e);

  // Compute 'kappa'
  norm = Mach*c*c+(1.-Mach)*vel*vel;
  norm *= _rho[_qp];
  Real kappa_e = std::max(res, jump);
  kappa_e *= h*h/norm;
  _kappa[_qp] = _is_first_order_visc || _t_step==1 ? _visc_max[_qp] :  std::min(_visc_max[_qp], kappa_e);

////    Real Mach2 = _areViscEqual ? Mach*Mach: 1.;
//    Real Mach2 = Mach*Mach;
//    Real fct_of_mach = Mach;
//    switch (_fct_of_mach_type) {
//        case MACH:
//            fct_of_mach = std::min(Mach, 1.);
//            break;
//        case SQRT_MACH:
//            fct_of_mach = std::min(std::sqrt(Mach), 1.);
//            break;
//        case FCT_OF_MACH:
//            fct_of_mach = std::min(Mach*std::sqrt(4+(1.-Mach*Mach)*(1.-Mach*Mach)) / (1.+Mach*Mach),1.);
//            break;
//        default:
//            mooseError("The function with name: \"" << _fct_of_mach_name << "\" is not supported in the \"ComputeViscosityCoefficient\" type of material.");
//    }
//    
//    // Postprocessors:
//    Real rhov2_pps = std::max(getPostprocessorValueByName(_rhov2_pps_name), eps);
//    Real alpha_var = getPostprocessorValueByName(_alpha_pps_name);
//    
//    // Initialyze some variables used in the switch statement:
//    Real weight0, weight1, weight2;
//    Real mu_e, kappa_e, beta_e;
//    Real jump, residual, norm, sigma;
//
//        case ENTROPY:
//            // Compute the weights for BDF2
//            if (_t_step > 1)
//            {
//                weight0 = (2.*_dt+_dt_old)/(_dt*(_dt+_dt_old));
//                weight1 = -(_dt+_dt_old)/(_dt*_dt_old);
//                weight2 = _dt/(_dt_old*(_dt+_dt_old));
//            }
//            else
//            {
//                weight0 =  1. / _dt;
//                weight1 = -1. / _dt;
//                weight2 = 0.;
//            }
//            
//        /** Compute viscosity coefficient for void fraction equation: **/
//            residual = 0.;
//            residual = _velI[_qp]*_grad_vf_liq[_qp];
//            residual += (weight0*_vf_liq[_qp]+weight1*_vf_liq_old[_qp]+weight2*_vf_liq_older[_qp]);
//            residual *= _Ce;
//            if (_is_jump_on)
//                jump = _Calpha*std::fabs(_velI[_qp].size()*_jump_grad_alpha[_qp]);
//            else
//                jump = _Calpha*std::fabs(_velI[_qp].size()*_grad_vf_liq[_qp](0));
////            norm = std::min(_vf_liq[_qp], std::fabs(1.-_vf_liq[_qp]));
//            norm = alpha_var;
//            beta_e = h*h*(std::fabs(residual)+jump) / norm;
////            beta_e += h*h*std::fabs(vel*_grad_area[_qp])/_area[_qp];
////            if (std::fabs(residual)>1e-3) {
////                std::cout<<"$$$$$$$$$$$$$$$$$"<<std::endl;
////                std::cout<<"alpha="<<_vf_liq[_qp]<<std::endl;
////                std::cout<<"alpha old="<<_vf_liq_old[_qp]<<std::endl;
////                std::cout<<"alpha older="<<_vf_liq_older[_qp]<<std::endl;
////                std::cout<<"grad="<<_grad_vf_liq[_qp](0)<<std::endl;
////                std::cout<<"vel="<<_velI[_qp]<<std::endl;
////                std::cout<<"residual="<<residual<<std::endl;
////                std::cout<<"pps="<<norm<<std::endl;
////            }
//        /** Compute viscosity coefficient for continuity, momentum and energy equations: **/
//            // Entropy residual:
//            residual = 0.;
//            residual = vel*_grad_press[_qp];
//            residual += (weight0*_pressure[_qp]+weight1*_pressure_old[_qp]+weight2*_pressure_older[_qp]);
//            residual -= c*c*vel*_grad_rho[_qp];
//            residual -= c*c*(weight0*_rho[_qp]+weight1*_rho_old[_qp]+weight2*_rho_older[_qp]);
//            residual *= _Ce;
//            
//            if (_isShock) // non-isentropic flow.
//            {
//                // Compute the jumps for mu_e:
//                if (_is_jump_on)
//                    jump = _Cjump*vel.size()*std::max( _jump_grad_press[_qp], Mach2*c*c*_jump_grad_dens[_qp] );
//                else
//                    jump = _Cjump*vel.size()*std::max( _grad_press[_qp].size(), Mach2*c*c*_grad_rho[_qp].size() );
//                
//                // Compute high-order viscosity coefficient mu_e:
//                norm = 0.5 * std::max( (1.-Mach)*rhov2_pps, _rho[_qp]*std::min(vel.size_sq(), c*c) );
//                mu_e = h*h*(std::fabs(residual) + jump) / norm;
//                mu_e += h*h*_pressure[_qp] * std::fabs(vel * _grad_area[_qp]) / ( _area[_qp] * norm );
//                                
//                // Compute the jumps for kappa_e:
//                if (_is_jump_on)
//                    jump = _Cjump*vel.size()*std::max( _jump_grad_press[_qp], c*c*_jump_grad_dens[_qp] );
//                else
//                    jump = _Cjump*vel.size()*std::max( _grad_press[_qp].size(), c*c*_grad_rho[_qp].size() );
//
//                // Compute high-order viscosity coefficient kappa_e:
////                norm = 0.5*( std::fabs(1.-Mach)*_rho[_qp]*c*c + Mach*_rho[_qp]*std::min(vel.size_sq(), c*c) );
//                if (!_areViscEqual)
//                    norm = 0.5*_rho[_qp]*c*c;
//                kappa_e = h*h*(std::fabs(residual) + jump) / norm;
//                kappa_e += h*h*_pressure[_qp] * std::fabs(vel * _grad_area[_qp]) / ( _area[_qp] * norm );
//            }
//            else // with function sigma.
//            {
//                // Compute the jumps:
//                if (_is_jump_on)
//                    jump = _Cjump*vel.size()*std::max( _jump_grad_press[_qp], c*c*_jump_grad_dens[_qp] );
//                else
//                    jump = _Cjump*vel.size()*std::max( _grad_press[_qp].size(), c*c*_grad_rho[_qp].size() );
//                
//                // Compute the function sigma:
//                sigma = 0.5 * std::tanh(_a_coeff*(Mach-_Mthres));
//                sigma += 0.5 * std::abs(std::tanh(_a_coeff*(Mach-_Mthres)));
//                
//                // Compute mu_e:
//                norm = (1.-sigma) * _rho[_qp] * c * c;
//                norm += sigma * _rho[_qp] * vel.size_sq();
//                norm *= 0.5;
//                mu_e = h*h*(std::fabs(residual) + jump) / norm;
//                
//                // Compute kappa_e:
//                norm = 0.5 * _rho[_qp] * c * c;
//                kappa_e = h*h*(std::fabs(residual) + jump) / norm;
//                
////                // Compute high-order viscosity coefficients:
////                norm = 0.5*( std::fabs(1.-Mach)*_rho[_qp]*c*c + Mach*_rho[_qp]*vel.size_sq() );
////                kappa_e = h*h*std::max(std::fabs(residual), jump) / norm;
//////                kappa_e += h*h*_pressure[_qp] * std::fabs(vel * _grad_area[_qp]) / ( _area[_qp] * norm );
////                kappa_e += h*h*std::fabs(vel*_grad_area[_qp])/_area[_qp];
////                mu_e = kappa_e;
//            }
//            
//        /** Compute the viscosity coefficients: **/
//            if (_t_step == 1)
//            {
//                _mu[_qp] = 2*_kappa_max[_qp];
//                _kappa[_qp] = 2*_kappa_max[_qp];
//                _beta[_qp] = 2*_beta_max[_qp];
//            }
//            else
//            {
//                _beta[_qp] = std::min(_beta_max[_qp], beta_e);
//                _kappa[_qp] = std::min( _kappa_max[_qp], kappa_e );
//                _mu[_qp] = std::min( _kappa_max[_qp], mu_e );
//            }
//            break;
//        default:
//            mooseError("The viscosity type entered in the input file is not implemented.");
//            break;
//    }
}
