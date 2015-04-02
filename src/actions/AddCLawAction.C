

#include "AddCLawAction.h"
#include "AddVariableAction.h"
#include "MooseApp.h"
#include "FEProblem.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/fe.h"

template<>
InputParameters validParams<AddCLawAction>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());
  InputParameters params = validParams<Action>();
  params.addParam<std::string>("type", " ");
  params.addParam<MooseEnum>("family", families, "Specifies the family of FE shape functions to use for this variable");
  params.addParam<MooseEnum>("order", orders,  "Specifies the order of the FE shape function to use for this variable (additional orders not listed are allowed)");
  params.addParam<std::vector<NonlinearVariableName> >("variables", "rho momentum_x momentum_y momentum_z rhoe" "非线性变量");
  return params;
}

AddCLawAction::AddCLawAction(const std::string & name, InputParameters params) :
    Action(name, params),
	_variables(getParam<std::vector<NonlinearVariableName> >("variables")),
	_type(getParam<std::string>("type"))
{
}

void AddCLawAction::act()
{
	if(_current_task == "add_variable")
		AddVariable();
	else if(_current_task == "add_kernel")
		AddKernel();
	else if(_current_task == "add_dg_kernel")
		AddDGKernel();
}

void AddCLawAction::AddVariable()
{
	Real scale_factor = isParamValid("scaling") ? getParam<Real>("scaling") : 1;
	FEType fe_type(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
             Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family")));
	for (int i = 0; i < _variables.size(); ++i)
	{
		_problem->addVariable(_variables[i], fe_type, scale_factor);
	}
}

void AddCLawAction::AddKernel()
{
	std::string time_kernel_name = "TimeDerivative";
	InputParameters params = _factory.getValidParams(time_kernel_name);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addKernel(time_kernel_name, _variables[i] + "_time", params);
	}

	std::string cell_kernel_name = "NSCellKernel";
	params = _factory.getValidParams(cell_kernel_name);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addKernel(cell_kernel_name, _variables[i] + "_space", params);
	}
}

void AddCLawAction::AddDGKernel()
{
	std::string face_kernel_name = "NSFaceKernel";
	InputParameters params = _factory.getValidParams(face_kernel_name);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addDGKernel(face_kernel_name, _variables[i] + "_dg", params);
	}
}

void AddCLawAction::AddAuxVariable()
{
}

void AddCLawAction::AddAuxKernel()
{
}
