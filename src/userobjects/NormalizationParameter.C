#include "NormalizationParameter.h"
#include "MooseError.h"

template<>
InputParameters validParams<NormalizationParameter>()
{
  InputParameters params = validParams<UserObject>();
  // Constants
  params.addRequiredParam<Real>("M_threshold", "Coefficient for function used in the normalization parameter");
  params.addRequiredParam<Real>("a", "Coefficient for function used in the normalization parameter");
  // Name of postprocessors
  params.addRequiredParam<PostprocessorName>("velocity_pps_name", "pps for velocity");
  // Function type
  params.addRequiredParam<std::string>("funct_type", "function used in the def. of the normalization of the visc. coeff.");
  return params;
}

NormalizationParameter::NormalizationParameter(const std::string & name, InputParameters parameters) :
    GeneralUserObject(name, parameters),
    // Constant
    _M_thres(getParam<Real>("M_threshold")),
    _a(getParam<Real>("a")),
    // Get the names of functions postprocessing velocity and pressure values
    _vel_pps(getPostprocessorValue("velocity_pps_name")),
    // Function type
    _fcnt_type("Mach_fnct Tanh_fnct Sin_fnct Shock_fnct", getParam<std::string>("funct_type"))
{
}

NormalizationParameter::~NormalizationParameter()
{
  // Destructor, empty
}

void
NormalizationParameter::destroy()
{
}

Real
NormalizationParameter::compute(Real c2, Real vel) const
{
  Real Mach = std::fabs(vel) / std::sqrt(c2);
  Real func, norm;
  switch (_fcnt_type)
  {
    case Mach_fnct:
      func = std::min(Mach, 1.);
      norm = 0.5 * std::fabs((1. - func) * c2 + func * std::min(c2, vel * vel));
      break;

    case Tanh_fnct:
      func = std::tanh(_a * (Mach - _M_thres));
      func += std::fabs(std::tanh(_a * (Mach - _M_thres)));
      func *= 0.5;
      norm = 0.5 * std::fabs((1. - func) * c2 + func * std::min(c2, vel * vel));
      break;

    case Sin_fnct:
      if (Mach < _M_thres - _a)
        func = 0.;
      else if (Mach > _M_thres + _a)
        func = 1.;
      else
      {
        func = 1. + (Mach - _M_thres) / _a;
        func += std::sin(libMesh::pi * (Mach - _M_thres) / _a) / libMesh::pi;
        func *= 0.5;
      }
      norm = 0.5 * std::fabs((1. - func) * c2 + func * std::min(c2, vel * vel));
      break;

    case Shock_fnct:
      norm = std::max(std::min(c2, vel * vel), _vel_pps * _vel_pps);
      break;

    default:
      mooseError("'" << this->name() << "' Undefined function type to compute the normalization parameter used in the computation of the visc. coeff.");
      break;
  }
  return norm;
}
