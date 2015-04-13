

#pragma once

#include "InputParameters.h"
#include "Assembly.h"
using std::vector;

class CLawProblem;
class NonlinearSystem;

class CLawInterface
{
public:
	CLawInterface(InputParameters &parameter);
	~CLawInterface(){}

	void convertionTerm(RealVectorValue *inviscous_term, Real *uh);
	void diffusionTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh);
	void sourceTerm();

	int equationIndex(const std::string& var_name);

	MooseVariable & getVariable(const std::string var_name);
	MooseVariable & getVariable(int eq);
protected:
	CLawProblem &_claw_problem;
	NonlinearSystem &_nl;
	THREAD_ID _tid;
	vector<VariableName> _variables;
	int _n_equations;
	int _var_order;
};
