
#include "CFDKernelsAction.h"
#include "MooseApp.h"
#include "FEProblem.h"

template<>
InputParameters validParams<CFDKernelsAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<std::string>("type", "CellKernel类型");
  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "非线性变量");
  return params;
}

CFDKernelsAction::CFDKernelsAction(const std::string & name, InputParameters params) :
    Action(name, params),
	_variables(getParam<std::vector<NonlinearVariableName> >("variables"))
{
}

void CFDKernelsAction::act()
{
	std::string time_kernel_name = "TimeDerivative";
	InputParameters params = _factory.getValidParams(time_kernel_name);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addKernel(time_kernel_name, _variables[i] + "_time", params);
	}

	std::string cell_kernel_name = getParam<std::string>("type");
	params = _factory.getValidParams(cell_kernel_name);
    _app.parser().extractParams(_name, params);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addKernel(cell_kernel_name, _variables[i] + "_space", params);
	}
}

