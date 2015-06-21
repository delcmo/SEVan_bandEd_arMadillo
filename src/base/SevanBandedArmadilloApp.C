#include "SevanBandedArmadilloApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

// kernels
#include "SbaVolumeFraction.h"
#include "SbaMass.h"
#include "SbaMomentum.h"
#include "SbaEnergy.h"
#include "SbaArtificialDissipation.h"

// auxkernels
#include "DensityAux.h"
#include "InternalEnergyAux.h"
#include "PressureAux.h"
#include "VelocityAux.h"
#include "VolumeFractionAux.h"
#include "EntropyVolumeFractionAux.h"

// initial conditions
#include "SbaICs.h"

// boundary conditions
#include "SbaDirichletBC.h"
#include "SbaOutletStaticPressureBC.h"
#include "SbaStagnationPandTBC.h"

// materials
#include "ComputeViscosityCoefficient.h"
#include "InterfacialRelaxationTransfer.h"

// equation of state
#include "EquationOfState.h"
#include "StiffenedGasEquationOfState.h"

// include userobjects
#include "JumpGradientInterface.h"
#include "NormalizationParameter.h"

// postprocessors
#include "TimeStepCFL.h"
#include "InfiniteNormFromAverageValue.h"

template<>
InputParameters validParams<SevanBandedArmadilloApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

SevanBandedArmadilloApp::SevanBandedArmadilloApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  SevanBandedArmadilloApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  SevanBandedArmadilloApp::associateSyntax(_syntax, _action_factory);
}

SevanBandedArmadilloApp::~SevanBandedArmadilloApp()
{
}

extern "C" void SevanBandedArmadilloApp__registerApps() { SevanBandedArmadilloApp::registerApps(); }
void
SevanBandedArmadilloApp::registerApps()
{
  registerApp(SevanBandedArmadilloApp);
}

void
SevanBandedArmadilloApp::registerObjects(Factory & factory)
{
  // kernels
  registerKernel(SbaMass);
  registerKernel(SbaMomentum);
  registerKernel(SbaEnergy);
  registerKernel(SbaVolumeFraction);
  registerKernel(SbaArtificialDissipation);

  // auxkernels
  registerAux(DensityAux);
  registerAux(InternalEnergyAux);
  registerAux(PressureAux);
  registerAux(VelocityAux);
  registerAux(VolumeFractionAux);
  registerAux(EntropyVolumeFractionAux);

  // initial conditions
  registerInitialCondition(SbaICs);

  // boundary conditions
  registerBoundaryCondition(SbaDirichletBC);
  registerBoundaryCondition(SbaOutletStaticPressureBC);
  registerBoundaryCondition(SbaStagnationPandTBC);

  // materials
  registerMaterial(ComputeViscosityCoefficient);
  registerMaterial(InterfacialRelaxationTransfer);

  // equation of state
  registerUserObject(EquationOfState);
  registerUserObject(StiffenedGasEquationOfState);

  // userobjects
  registerUserObject(JumpGradientInterface);
  registerUserObject(NormalizationParameter);

  // postprocessors
  registerPostprocessor(TimeStepCFL);
  registerPostprocessor(InfiniteNormFromAverageValue);
}

void
SevanBandedArmadilloApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
