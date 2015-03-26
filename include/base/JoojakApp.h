#ifndef JOOJAKAPP_H
#define JOOJAKAPP_H

#include "MooseApp.h"

class JoojakApp;

template<>
InputParameters validParams<JoojakApp>();

class JoojakApp : public MooseApp
{
public:
  JoojakApp(const std::string & name, InputParameters parameters);
  virtual ~JoojakApp();

  void printHeader();
  virtual void run();
  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
  static void registerSAObjects(Factory & factory);
};

#endif /* JOOJAKAPP_H */
