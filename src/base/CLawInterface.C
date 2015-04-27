
#include "CLawInterface.h"
#include "CLawProblem.h"

CLawInterface::CLawInterface(InputParameters& parameter):
	_claw_problem(static_cast<CLawProblem&>(*parameter.get<FEProblem *>("_fe_problem"))),
	_nl(_claw_problem.getNonlinearSystem()),
	_tid(parameter.get<THREAD_ID>("_tid")),
	_variables(_nl.getVariableNames()),
	_n_equations(_variables.size()),
	_var_order(_claw_problem.getVariable(_tid, _variables[0]).order())
{
	if(_n_equations > 10)
		mooseError("当前支持方程个数不超过10");

	if(this->_n_equations != _claw_problem._n_equations)
		mooseError("参数文件中指定的变量个数和问题的变量数数目不同，检查Problem block下的参数" << _claw_problem.name());
}

//void CLawInterface::inviscousTerm(RealVectorValue *inviscous_term, Real *uh)
//{
//	_claw_problem.inviscousTerm(inviscous_term, uh);
//}
//
//void CLawInterface::viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh)
//{
//	_claw_problem.viscousTerm(viscous_term, uh, duh);
//}
//
//void CLawInterface::sourceTerm()
//{
////	_claw_base.sourceTerm();
//}

int CLawInterface::equationIndex(const std::string& var_name)
{
	return _claw_problem.equationIndex(var_name);
}

MooseVariable& CLawInterface::getVariable(const std::string var_name)
{
	return _claw_problem.getVariable(_tid, var_name);
}

MooseVariable& CLawInterface::getVariable(int eq)
{
	return _claw_problem.getVariable(_tid, _variables[eq]);
}
