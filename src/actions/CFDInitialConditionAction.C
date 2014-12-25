
#include "CFDInitialConditionAction.h"
#include "MooseApp.h"
#include "FEProblem.h"

template<>
InputParameters validParams<CFDInitialConditionAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<std::string>("type", "InitialCondition类型");
  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "非线性变量");
  return params;
}

CFDInitialConditionAction::CFDInitialConditionAction(const std::string & name, InputParameters params) :
    Action(name, params),
	_variables(getParam<std::vector<NonlinearVariableName> >("variables"))
{
}

void CFDInitialConditionAction::act()
{
	std::string init_cond_name = getParam<std::string>("type");
    InputParameters params = _factory.getValidParams(init_cond_name);
    _app.parser().extractParams(_name, params);
	for (int i = 0; i < _variables.size(); ++i)
	{
	    params.set<VariableName>("variable") = _variables[i];
	    _problem->addInitialCondition(init_cond_name, _variables[i]+"_ic", params);
	}
}

