
#pragma once

#include "Action.h"

class AddCLawAction : public Action
{
public:
  AddCLawAction(const std::string & name, InputParameters params);
  virtual void act();

private:
  void AddVariable();
  void AddKernel();
  void AddDGKernel();
  void AddAuxVariable();
  void AddAuxKernel();
protected:
  std::vector<NonlinearVariableName> _variables;
  std::string _type;

};

template<>
InputParameters validParams<AddCLawAction>();
