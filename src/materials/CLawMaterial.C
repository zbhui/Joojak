
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
		_var_order(_claw_problem.getVariable(_tid, _variables[0]).order()),
		_current_elem_volume(_assembly.elemVolume()),
		_neighbor_elem_volume(_assembly.neighborVolume()),
		_current_side_volume(_assembly.sideElemVolume())
{
	_aux_variables = _claw_problem._aux_variables;
	for(int i = 0; i < _aux_variables.size(); ++i)
		_variables.push_back(_aux_variables[i]);
	_n_variables = _variables.size();

	for (int ivar = 0; ivar < _n_variables; ++ivar)
	{
		MooseVariable &val = _claw_problem.getVariable(_tid, _variables[ivar]);
		_uh.push_back(_is_implicit ? &val.sln() : &val.slnOld());
		_uh_neighbor.push_back(_is_implicit ? &val.slnNeighbor() : &val.slnOldNeighbor());
		_grad_uh.push_back(_is_implicit ? &val.gradSln(): &val.gradSlnOld());
		_grad_uh_neighbor.push_back(_is_implicit ? &val.gradSlnNeighbor(): &val.gradSlnOldNeighbor());

		addMooseVariableDependency(&_fe_problem.getVariable(_tid, _variables[ivar]));
	}
}
