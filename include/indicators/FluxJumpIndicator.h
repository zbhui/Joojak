
#pragma once

#include "InternalSideIndicator.h"
#include "CLawFaceMaterial.h"

class CLawProblem;
class CFDProblem;

class FluxJumpIndicator :
public InternalSideIndicator
{
public:
  FluxJumpIndicator(const std::string & name, InputParameters parameters);
  virtual ~FluxJumpIndicator(){};

protected:
	CLawProblem &_claw_problem;
	CFDProblem &_cfd_problem;
	NonlinearSystem &_nl;
	THREAD_ID _tid;
	vector<VariableName> _variables;
	vector<VariableName> _aux_variables;
	int _n_variables;
	int _var_order;

	vector<VariableValue*> _uh;
	vector<VariableValue*> _uh_neighbor;
	vector<VariableGradient*> _grad_uh;
	vector<VariableGradient*> _grad_uh_neighbor;

	bool _is_implicit;
  virtual Real computeQpIntegral();
  void computeIndicator();
  void finalize();
};

template<>
InputParameters validParams<FluxJumpIndicator>();
