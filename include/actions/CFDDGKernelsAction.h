
#pragma once

#include "Action.h"

class CFDDGKernelsAction;

template<>
InputParameters validParams<CFDDGKernelsAction>();

class CFDDGKernelsAction : public Action
{
public:
  CFDDGKernelsAction(const std::string & name, InputParameters params);
  virtual void act();

protected:
  std::vector<NonlinearVariableName> _variables;
};
