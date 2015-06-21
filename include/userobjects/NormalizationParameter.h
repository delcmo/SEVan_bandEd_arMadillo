#ifndef NORMALIZATIONPARAMETER_H
#define NORMALIZATIONPARAMETER_H

#include "GeneralUserObject.h"

// Forward Declarations
class NormalizationParameter;

template<>
InputParameters validParams<NormalizationParameter>();

class NormalizationParameter : public GeneralUserObject
{
public:
  // Constructor
  NormalizationParameter(const std::string & name, InputParameters parameters);

  // Destructor
  virtual ~NormalizationParameter();

  /**
   * Called when this object needs to compute something.
   */
  virtual void execute() {}

  /**
   * Called before execute() is ever called so that data can be cleared.
   */
  virtual void initialize(){}


  virtual void destroy();

  virtual void finalize() {};

  // The interface for derived NormalizationParameter objects to implement...
  virtual Real compute(Real c2 = 0., Real vel = 0.) const;

protected:
  // Constants
  Real _M_thres;
  Real _a;

  // Pps for velocity
  const PostprocessorValue & _vel_pps;

  // Function type
  enum FnctType
  {
    Mach_fnct = 0,
    Tanh_fnct = 1,
    Sin_fnct = 2,
    Shock_fnct = 3
  };
  MooseEnum _fcnt_type;
};


#endif // NORMALIZATIONPARAMETER_H
