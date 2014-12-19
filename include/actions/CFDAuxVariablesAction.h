
#pragma once

#include "Action.h"

class CFDAuxVariablesAction;

template<>
InputParameters validParams<CFDAuxVariablesAction>();

class CFDAuxVariablesAction : public Action
{
public:
	CFDAuxVariablesAction(const std::string & name, InputParameters params);
  virtual void act();

protected:
  std::vector<NonlinearVariableName> _variables;
  std::vector<AuxVariableName> _aux_variables;
};
