#include "SevanBandedArmadilloApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

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
}

void
SevanBandedArmadilloApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
