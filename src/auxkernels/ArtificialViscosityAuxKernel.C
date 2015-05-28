
#include "ArtificialViscosityAuxKernel.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<ArtificialViscosityAuxKernel>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("indicator", "indicator");
  params.addRequiredCoupledVar("marker", "marker");
  return params;
}

ArtificialViscosityAuxKernel::ArtificialViscosityAuxKernel(const std::string & name, InputParameters parameter) :
	AuxKernel(name, parameter),
	_cfd_problem(static_cast<CFDProblem&>(*parameter.get<FEProblem *>("_fe_problem"))),
	_nl(_cfd_problem.getNonlinearSystem()),
	_tid(parameter.get<THREAD_ID>("_tid")),
	_variables(_nl.getVariableNames()),
	_n_equations(_variables.size()),
	_var_order(_cfd_problem.getVariable(_tid, _variables[0]).order()),
	_indicator(coupledValue("indicator"))
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		MooseVariable &val = _cfd_problem.getVariable(_tid, _variables[eq]);
		_uh.push_back(&val.sln());
	}
}

Real ArtificialViscosityAuxKernel::computeValue()
{
	Real uh[10];
	std::string var_name = _var.name();
	valueAtCellPoint(uh);

	return _indicator[_qp]*_current_elem->hmax();
//	return _cfd_problem.computeAuxValue(var_name, uh);

}

void ArtificialViscosityAuxKernel::valueAtCellPoint(Real *uh)
{
	for (size_t eq = 0; eq < _n_equations; ++eq)
	{
		uh[eq] = (*_uh[eq])[_qp];
	}
}
