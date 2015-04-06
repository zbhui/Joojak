
#pragma once

#include "Action.h"

class CFDInitialConditionAction;

template<>
InputParameters validParams<CFDInitialConditionAction>();

class CFDInitialConditionAction : public Action
{
public:
  CFDInitialConditionAction(const std::string & name, InputParameters params);
  virtual void act();

protected:
};
