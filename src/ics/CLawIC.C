
#include "CLawIC.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawIC>()
{
  InputParameters params = validParams<MultiInitialCondition>();
  return params;
}

CLawIC::CLawIC(const std::string & name, InputParameters parameters) :
		MultiInitialCondition(name, parameters),
		_claw_problem(static_cast<CLawProblem&>(_fe_problem))
{
}

Real CLawIC::value(int component, const Point & p)
{
	return _claw_problem.initialCondition(p, component);
}
