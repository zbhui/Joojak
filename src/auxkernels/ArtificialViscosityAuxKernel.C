
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
		_grad_uh.push_back(_is_implicit ? &val.gradSln(): &val.gradSlnOld());

	}
}

Real ArtificialViscosityAuxKernel::computeValue()
{
	Real uh[10];
	RealVectorValue duh[10];
	std::string var_name = _var.name();
	valueAtCellPoint(uh, duh);

	Real viscosity = _indicator[_qp]*_current_elem->hmax();
//	if (_indicator[_qp] > 0)
//	  std :: cout << duh[4].size()/uh[4] << std::endl;

//	viscosity *= (duh[4].size()/uh[4] * _current_elem->hmax());

//	if(_t_step < 100)
//		return 500*viscosity;
//	else
//		return 10*viscosity;

	return viscosity;

}

void ArtificialViscosityAuxKernel::valueAtCellPoint(Real *uh, RealVectorValue *duh)
{
	for (size_t eq = 0; eq < _n_equations; ++eq)
	{
		uh[eq] = (*_uh[eq])[_qp];
		duh[eq] = (*_grad_uh[eq])[_qp];
	}
}
