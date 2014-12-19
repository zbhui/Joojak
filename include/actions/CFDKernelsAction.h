
#pragma once

#include "Action.h"

class CFDKernelsAction;

template<>
InputParameters validParams<CFDKernelsAction>();

class CFDKernelsAction : public Action
{
public:
  CFDKernelsAction(const std::string & name, InputParameters params);
  virtual void act();

protected:
  std::vector<NonlinearVariableName> _variables;

  void addKernels();
};
