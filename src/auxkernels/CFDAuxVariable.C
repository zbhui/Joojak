
#include "CFDAuxVariable.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<CFDAuxVariable>()
{
  InputParameters params = validParams<AuxKernel>();
  return params;
}

CFDAuxVariable::CFDAuxVariable(const std::string & name, InputParameters parameter) :
    AuxKernel(name, parameter),
	_cfd_problem(static_cast<CFDProblem&>(*parameter.get<FEProblem *>("_fe_problem"))),
	_nl(_cfd_problem.getNonlinearSystem()),
	_tid(parameter.get<THREAD_ID>("_tid")),
	_variables(_nl.getVariableNames()),
	_n_equations(_variables.size()),
	_var_order(_cfd_problem.getVariable(_tid, _variables[0]).order())
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		MooseVariable &val = _cfd_problem.getVariable(_tid, _variables[eq]);
		_uh.push_back(&val.sln());
	}
}

Real CFDAuxVariable::computeValue()
{
	Real uh[10];
	std::string var_name = _var.name();
	valueAtCellPoint(uh);

	return _cfd_problem.computeAuxValue(var_name, uh);

}

void CFDAuxVariable::valueAtCellPoint(Real *uh)
{
	for (size_t eq = 0; eq < _n_equations; ++eq)
	{
		uh[eq] = (*_uh[eq])[_qp];
	}
}
