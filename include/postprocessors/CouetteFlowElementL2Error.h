
#pragma once

#include "ElementIntegralPostprocessor.h"

class CouetteFlowProblem;

using std::vector;

class CouetteFlowElementL2Error :
public ElementIntegralPostprocessor
{
public:
	CouetteFlowElementL2Error(const std::string & name, InputParameters parameters);

	virtual Real getValue();

protected:
	virtual Real computeQpIntegral();

private:
	CouetteFlowProblem & _couette_problem;
	NonlinearSystem &_nl;
	THREAD_ID _tid;
	vector<VariableName> _variables;
    int _n_equations;

	std::vector<VariableValue*> _uh;


};

template<>
InputParameters validParams<CouetteFlowElementL2Error>();
