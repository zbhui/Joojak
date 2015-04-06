
#include "CFDInitialConditionAction.h"
#include "MooseApp.h"
#include "FEProblem.h"

template<>
InputParameters validParams<CFDInitialConditionAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<std::string>("type", "InitialCondition类型");
  return params;
}

CFDInitialConditionAction::CFDInitialConditionAction(const std::string & name, InputParameters params) :
    Action(name, params)
{
}

void CFDInitialConditionAction::act()
{
	std::vector<VariableName> var = _problem->getNonlinearSystem().getVariableNames();
	std::string init_cond_name = getParam<std::string>("type");
    InputParameters params = _factory.getValidParams(init_cond_name);
    _app.parser().extractParams(_name, params);
	for (int i = 0; i < var.size(); ++i)
	{
	    params.set<VariableName>("variable") = var[i];
	    _problem->addInitialCondition(init_cond_name, var[i]+"_ic", params);
	}
}

