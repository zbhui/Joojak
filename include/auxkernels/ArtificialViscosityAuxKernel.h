
#pragma once

#include "AuxKernel.h"

class CFDProblem;
using std::vector;

class ArtificialViscosityAuxKernel :
public AuxKernel
{
public:
	ArtificialViscosityAuxKernel(const std::string & name, InputParameters parameters);

protected:

  CFDProblem &_cfd_problem;
  NonlinearSystem &_nl;
  THREAD_ID _tid;
  vector<VariableName> _variables;
  int _n_equations;
  int _var_order;

  vector<VariableValue*> _uh;
  vector<VariableGradient*> _grad_uh;
  VariableValue &_indicator;
  virtual Real computeValue();
  void valueAtCellPoint(Real *uh, RealVectorValue *duh);
};

template<>
InputParameters validParams<ArtificialViscosityAuxKernel>();
