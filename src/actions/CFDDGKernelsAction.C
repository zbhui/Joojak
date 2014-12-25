
#include "CFDDGKernelsAction.h"
#include "MooseApp.h"
#include "FEProblem.h"

template<>
InputParameters validParams<CFDDGKernelsAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<std::string>("type", "FaceKernel类型");
  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "非线性变量");
  return params;
}

CFDDGKernelsAction::CFDDGKernelsAction(const std::string & name, InputParameters params) :
    Action(name, params),
	_variables(getParam<std::vector<NonlinearVariableName> >("variables"))
{
}

void CFDDGKernelsAction::act()
{
	std::string face_kernel_name = getParam<std::string>("type");
	InputParameters params = _factory.getValidParams(face_kernel_name);
    _app.parser().extractParams(_name, params);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addDGKernel(face_kernel_name, _variables[i] + "_dg", params);
	}
}

