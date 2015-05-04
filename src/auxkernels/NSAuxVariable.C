
#include "NSAuxVariable.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<NSAuxVariable>()
{
  InputParameters params = validParams<AuxKernel>();
  return params;
}

NSAuxVariable::NSAuxVariable(const std::string & name, InputParameters parameter) :
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

Real NSAuxVariable::computeValue()
{
	Real uh[10];
	std::string var_name = _var.name();
	valueAtCellPoint(uh);

	if(var_name == "pressure")
		return _cfd_problem.pressure(uh);
	if(var_name == "mach")
		return _cfd_problem.localMach(uh);
	if(var_name == "velocity_x")
		return uh[1]/uh[0];
	if(var_name == "velocity_y")
		return uh[2]/uh[0];
	if(var_name == "velocity_z")
		return uh[3]/uh[0];

	mooseError(var_name << "辅助变量名不存在");
	return 0.;

}

void NSAuxVariable::valueAtCellPoint(Real *uh)
{
	for (size_t eq = 0; eq < _n_equations; ++eq)
	{
		uh[eq] = (*_uh[eq])[_qp];
	}
}
