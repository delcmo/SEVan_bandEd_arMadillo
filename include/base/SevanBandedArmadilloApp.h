#ifndef SEVAN_BANDED_ARMADILLOAPP_H
#define SEVAN_BANDED_ARMADILLOAPP_H

#include "MooseApp.h"

class SevanBandedArmadilloApp;

template<>
InputParameters validParams<SevanBandedArmadilloApp>();

class SevanBandedArmadilloApp : public MooseApp
{
public:
  SevanBandedArmadilloApp(const std::string & name, InputParameters parameters);
  virtual ~SevanBandedArmadilloApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* SEVAN_BANDED_ARMADILLOAPP_H */
