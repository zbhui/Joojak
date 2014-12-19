
#include "CFDAuxVariablesAction.h"
#include "AddVariableAction.h"
#include "MooseApp.h"
#include "FEProblem.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/fe.h"

template<>
InputParameters validParams<CFDAuxVariablesAction>()
{
	MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
	MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());

	InputParameters params = validParams<Action>();
	params.addRequiredParam<std::string>("type", "AuxKernel类型");
	params.addParam<MooseEnum>("family", families, "Specifies the family of FE shape functions to use for this variable");
	params.addParam<MooseEnum>("order", orders,  "Specifies the order of the FE shape function to use for this variable (additional orders not listed are allowed)");
	params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "非线性变量");
	params.addParam<std::vector<AuxVariableName> >("aux_variables", "辅助变量名");
	return params;
}

CFDAuxVariablesAction::CFDAuxVariablesAction(const std::string & name, InputParameters params) :
    		Action(name, params),
			_variables(getParam<std::vector<NonlinearVariableName> >("variables")),
			_aux_variables(getParam<std::vector<AuxVariableName> >("aux_variables"))
{
}

void CFDAuxVariablesAction::act()
{
	if(_current_task == "add_aux_variable")
	{
		FEType fe_type(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
				Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family")));

		for (int i = 0; i < _aux_variables.size(); ++i)
		{
			_problem->addAuxVariable(_aux_variables[i], fe_type);
		}
	}

	else if(_current_task == "add_aux_kernel")
	{
		std::vector<VariableName> var_name;
		for (int i = 0; i < _variables.size(); ++i)
			var_name.push_back(_variables[i]);

		const std::string aux_kernel_name = getParam<std::string>("type");
		InputParameters params = _factory.getValidParams(aux_kernel_name);
		params.set<std::vector<VariableName> >("variables") = var_name;
		_app.parser().extractParams(_name, params);

		for (int i = 0; i < _aux_variables.size(); ++i)
		{
			params.set<AuxVariableName>("variable") = _aux_variables[i];
			_problem->addAuxKernel(aux_kernel_name, _aux_variables[i], params);
		}
	}
	else
	{
		mooseError("unknown task.");
	}
}

