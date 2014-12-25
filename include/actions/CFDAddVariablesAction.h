
#pragma once

#include "Action.h"

class CFDAddVariablesAction;

template<>
InputParameters validParams<CFDAddVariablesAction>();

class CFDAddVariablesAction : public Action
{
public:
  CFDAddVariablesAction(const std::string & name, InputParameters params);
  virtual void act();

protected:
  std::vector<NonlinearVariableName> _variables;
};
