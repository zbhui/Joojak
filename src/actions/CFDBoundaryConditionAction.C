
#include "CFDBoundaryConditionAction.h"
#include "MooseApp.h"
#include "FEProblem.h"

template<>
InputParameters validParams<CFDBoundaryConditionAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<std::string>("type", "InitialCondition类型");
  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "非线性变量");
  params.addParam<std::vector<BoundaryName> >("boundary", "CFD边界");
  return params;
}

CFDBoundaryConditionAction::CFDBoundaryConditionAction(const std::string & name, InputParameters params) :
    Action(name, params),
	_variables(getParam<std::vector<NonlinearVariableName> >("variables"))
{
}

void CFDBoundaryConditionAction::act()
{
	std::string boun_cond_name = getParam<std::string>("type");
    InputParameters params = _factory.getValidParams(boun_cond_name);
    params.set<std::vector<BoundaryName> >("boundary") = getParam<std::vector<BoundaryName> >("boundary");
    _app.parser().extractParams(_name, params);
	for (int i = 0; i < _variables.size(); ++i)
	{
	    params.set<NonlinearVariableName>("variable") = _variables[i];
	    _problem->addBoundaryCondition(boun_cond_name, _variables[i]+"_bc", params);
	}
}

