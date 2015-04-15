
#pragma once

#include "Action.h"

class AddCLawAction : public Action
{
public:
  AddCLawAction(const std::string & name, InputParameters params);
  virtual void act();

private:
  void addVariable();
  void addKernel();
  void addDGKernel();
  void addAuxVariable();
  void addAuxKernel();
  void addBoundaryCondition();
protected:
  std::vector<NonlinearVariableName> _variables;
};

template<>
InputParameters validParams<AddCLawAction>();
