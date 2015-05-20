#include "IsoVortexElementL2Error.h"
#include "IsoVortexProblem.h"

template<>
InputParameters validParams<IsoVortexElementL2Error>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  return params;
}

IsoVortexElementL2Error::IsoVortexElementL2Error(const std::string & name, InputParameters parameters) :
	ElementIntegralPostprocessor(name, parameters),
	_isovortex_problem(static_cast<IsoVortexProblem&>(_fe_problem)),
	_nl(_isovortex_problem.getNonlinearSystem()),
	_tid(parameters.get<THREAD_ID>("_tid")),
	_variables(_nl.getVariableNames()),
	_n_equations(_variables.size())
{
	for (int eq = 0; eq < _isovortex_problem._n_equations; ++eq)
	{
		MooseVariable &val = _isovortex_problem.getVariable(_tid, _variables[eq]);
		_uh.push_back(_is_implicit ? &val.sln() : &val.slnOld());
	}
}

Real IsoVortexElementL2Error::getValue()
{
  return std::sqrt(ElementIntegralPostprocessor::getValue());
}

Real IsoVortexElementL2Error::computeQpIntegral()
{
	Real uh[5];
	for (int eq = 0; eq < _isovortex_problem._n_equations; ++eq)
		uh[eq] = (*_uh[eq])[_qp];

	Real err = uh[0] - _isovortex_problem.valueExact(_t, _q_point[_qp] , 0);
	return err*err;
}
