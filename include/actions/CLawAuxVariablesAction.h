
#pragma once

#include "Action.h"

class CLawAuxVariablesAction : public Action
{
public:
  CLawAuxVariablesAction(const std::string & name, InputParameters params);
  virtual void act();

protected:
  std::vector<AuxVariableName> _aux_variables;
};

template<>
InputParameters validParams<CLawAuxVariablesAction>();
