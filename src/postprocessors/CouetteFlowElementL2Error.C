#include "CouetteFlowElementL2Error.h"
#include "CouetteFlowProblem.h"

template<>
InputParameters validParams<CouetteFlowElementL2Error>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  return params;
}

CouetteFlowElementL2Error::CouetteFlowElementL2Error(const std::string & name, InputParameters parameters) :
	ElementIntegralPostprocessor(name, parameters),
	_couette_problem(static_cast<CouetteFlowProblem&>(_fe_problem)),
	_nl(_couette_problem.getNonlinearSystem()),
	_tid(parameters.get<THREAD_ID>("_tid")),
	_variables(_nl.getVariableNames()),
	_n_equations(_variables.size())
{
	for (int eq = 0; eq < _couette_problem._n_equations; ++eq)
	{
		MooseVariable &val = _couette_problem.getVariable(_tid, _variables[eq]);
		_uh.push_back(_is_implicit ? &val.sln() : &val.slnOld());
	}
}

Real CouetteFlowElementL2Error::getValue()
{
  return std::sqrt(ElementIntegralPostprocessor::getValue());
}

Real CouetteFlowElementL2Error::computeQpIntegral()
{
	Real uh[5];
	for (int eq = 0; eq < _couette_problem._n_equations; ++eq)
		uh[eq] = (*_uh[eq])[_qp];

	Real err = uh[0] - _couette_problem.valueExact(_t, _q_point[_qp] , 0);
	return err*err;
}
