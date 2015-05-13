
#include "CLawMaterial.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawMaterial>()
{
  InputParameters params = validParams<Material>();
  return params;
}

CLawMaterial::CLawMaterial(const std::string & name, InputParameters parameter):
		Material(name, parameter),
		_claw_problem(static_cast<CLawProblem&>(_fe_problem)),
		_nl(_claw_problem.getNonlinearSystem()),
		_tid(parameter.get<THREAD_ID>("_tid")),
		_variables(_nl.getVariableNames()),
		_n_equations(_variables.size()),
		_n_variables(coupledComponents("variables")),
		_var_order(_claw_problem.getVariable(_tid, _variables[0]).order())
{
}
