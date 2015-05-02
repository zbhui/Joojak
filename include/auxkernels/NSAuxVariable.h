
#pragma once

#include "AuxKernel.h"
#include "CLawInterface.h"

class CFDProblem;

class NSAuxVariable :
public AuxKernel,
public CLawInterface
{
public:
  NSAuxVariable(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
  void valueAtCellPoint(Real *uh);

  vector<VariableValue*> _uh;
  CFDProblem &_cfd_problem;
};

template<>
InputParameters validParams<NSAuxVariable>();
