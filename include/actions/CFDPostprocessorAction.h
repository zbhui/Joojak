
#pragma once

#include "Action.h"

class CFDPostprocessorAction;

template<>
InputParameters validParams<CFDPostprocessorAction>();

class CFDPostprocessorAction : public Action
{
public:
  CFDPostprocessorAction(const std::string & name, InputParameters params);
  virtual void act();
};
