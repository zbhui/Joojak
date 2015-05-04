
#pragma once

#include "AuxKernel.h"
//#include "CLawInterface.h"

class CFDProblem;
using std::vector;

class NSAuxVariable :
public AuxKernel
{
public:
  NSAuxVariable(const std::string & name, InputParameters parameters);

protected:

  CFDProblem &_cfd_problem;
  NonlinearSystem &_nl;
  THREAD_ID _tid;
  vector<VariableName> _variables;
  int _n_equations;
  int _var_order;

  vector<VariableValue*> _uh;

  virtual Real computeValue();
  void valueAtCellPoint(Real *uh);
};

template<>
InputParameters validParams<NSAuxVariable>();
