
#include "CFDAddVariablesAction.h"
#include "AddVariableAction.h"
#include "MooseApp.h"
#include "FEProblem.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/fe.h"

template<>
InputParameters validParams<CFDAddVariablesAction>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());
  InputParameters params = validParams<Action>();
  params.addParam<MooseEnum>("family", families, "Specifies the family of FE shape functions to use for this variable");
  params.addParam<MooseEnum>("order", orders,  "Specifies the order of the FE shape function to use for this variable (additional orders not listed are allowed)");
  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "非线性变量");
  return params;
}

CFDAddVariablesAction::CFDAddVariablesAction(const std::string & name, InputParameters params) :
    Action(name, params),
	_variables(getParam<std::vector<NonlinearVariableName> >("variables"))
{
}

void CFDAddVariablesAction::act()
{
	Real scale_factor = isParamValid("scaling") ? getParam<Real>("scaling") : 1;
	FEType fe_type(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
             Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family")));
	for (int i = 0; i < _variables.size(); ++i)
	{
		_problem->addVariable(_variables[i], fe_type, scale_factor);
	}

}
