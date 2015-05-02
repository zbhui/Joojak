
#include "NSAuxVariable.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<NSAuxVariable>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredCoupledVar("variables", "守恒变量");
  return params;
}

NSAuxVariable::NSAuxVariable(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
	CLawInterface(parameters),
	_cfd_problem(static_cast<CFDProblem&>(_claw_problem))
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		MooseVariable &val = getVariable(eq);
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
	size_t n_equation = coupledComponents("variables");
	for (size_t eq = 0; eq < n_equation; ++eq)
	{
//		uh[eq] = coupledValue("variables", eq)[_qp];
		uh[eq] = (*_uh[eq])[_qp];
	}
}
