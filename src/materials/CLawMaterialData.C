
#include "CLawMaterialData.h"
#include "CLawProblem.h"

CLawMaterialData::CLawMaterialData(const std::string & name, InputParameters parameters) :
    _subproblem(*parameters.get<SubProblem *>("_subproblem")),
    _fe_problem(*parameters.get<FEProblem *>("_fe_problem")),
	_claw_problem(static_cast<CLawProblem&>(_fe_problem)),
    _assembly(_subproblem.assembly(_tid)),
	_nl(_claw_problem.getNonlinearSystem()),
	_var_name(_nl.getVariableNames()),
	_n_variables(_var_name.size()),
    _tid(parameters.get<THREAD_ID>("_tid")),
    _bnd(parameters.get<bool>("_bnd")),
    _neighbor(parameters.get<bool>("_neighbor"))
{
	for (int eq = 0; eq < _n_variables; ++eq)
	{
		_moose_var.push_back(&_claw_problem.getVariable(_tid, _var_name[eq]));
	}
}

CLawMaterialData::~CLawMaterialData()
{
}

void CLawMaterialData::computeMaterial()
{
//	_claw_problem.computeMaterial(*this);
}

bool CLawMaterialData::isBoundary()
{
	return (_bnd );
}

bool CLawMaterialData::isElement()
{
	return !(_bnd || _neighbor);
}

bool CLawMaterialData::isSideSet()
{
	return (_bnd && _neighbor);
}
