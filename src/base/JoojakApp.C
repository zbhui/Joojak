#include "JoojakApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

template<>
InputParameters validParams<JoojakApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

JoojakApp::JoojakApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  JoojakApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  JoojakApp::associateSyntax(_syntax, _action_factory);
}

JoojakApp::~JoojakApp()
{
}

void
JoojakApp::registerApps()
{
  registerApp(JoojakApp);
}

void
JoojakApp::registerObjects(Factory & factory)
{
}

void
JoojakApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
