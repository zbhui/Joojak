
#pragma once

#include "ElementIntegralPostprocessor.h"

class IsoVortexProblem;

using std::vector;

class IsoVortexElementL2Error :
public ElementIntegralPostprocessor
{
public:
	IsoVortexElementL2Error(const std::string & name, InputParameters parameters);

	virtual Real getValue();

protected:
	virtual Real computeQpIntegral();

private:
	IsoVortexProblem & _isovortex_problem;
	NonlinearSystem &_nl;
	THREAD_ID _tid;
	vector<VariableName> _variables;
    int _n_equations;

	std::vector<VariableValue*> _uh;


};

template<>
InputParameters validParams<IsoVortexElementL2Error>();
