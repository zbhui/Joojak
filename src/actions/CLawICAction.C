
#include "CLawICAction.h"

#include "MooseApp.h"
#include "FEProblem.h"

template<>
InputParameters validParams<CLawICAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<std::string>("type", "InitialCondition类型");
  return params;
}

CLawICAction::CLawICAction(const std::string & name, InputParameters params) :
    Action(name, params)
{
}

void CLawICAction::act()
{
	std::vector<VariableName> var = _problem->getNonlinearSystem().getVariableNames();
	std::string init_cond_name = getParam<std::string>("type");
    InputParameters params = _factory.getValidParams(init_cond_name);
    _app.parser().extractParams(_name, params);
	for (int i = 0; i < var.size(); ++i)
	{
	    params.set<VariableName>("variable") = var[i];
	    params.set<int>("component") = i;
	    _problem->addInitialCondition(init_cond_name, var[i]+"_ic", params);
	}
}

