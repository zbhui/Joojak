
#pragma once

#include "Action.h"

class CFDBoundaryConditionAction;

template<>
InputParameters validParams<CFDBoundaryConditionAction>();

class CFDBoundaryConditionAction : public Action
{
public:
  CFDBoundaryConditionAction(const std::string & name, InputParameters params);
  virtual void act();

protected:
  std::vector<NonlinearVariableName> _variables;
};
