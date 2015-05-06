
#include "AddMultiAuxVariableAction.h"
#include "AddKernelAction.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "EigenSystem.h"

#include "libmesh/libmesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"

template<>
InputParameters validParams<AddMultiAuxVariableAction>()
{
  InputParameters params = validParams<AddAuxVariableAction>();
  params += validParams<AddKernelAction>();
  params.addRequiredParam<std::vector<AuxVariableName> >("variables", "非线性变量");
  params.addRequiredParam<std::string>("type", "A string representing the Moose Object that will be built by this Action");

  return params;
}

AddMultiAuxVariableAction::AddMultiAuxVariableAction(const std::string & name, InputParameters params) :
	AddAuxVariableAction(name, params),
	_type(getParam<std::string>("type")),
	_variables(getParam<std::vector<AuxVariableName> >("variables"))
{
}

void AddMultiAuxVariableAction::act()
{
	for (int i = 0; i < _variables.size(); ++i)
	{
		std::string var_name = _variables[i];
		if (_current_task == "add_aux_variable")
			addAuxVariable(var_name);
		else if(_current_task == "add_aux_kernel")
			addAuxKernel(var_name);
	}
}

void AddMultiAuxVariableAction::addAuxKernel(std::string var_name)
{
	InputParameters params = _factory.getValidParams(_type);
	_app.parser().extractParams(_name, params);
	params.set<AuxVariableName>("variable") = var_name;
	_problem->addAuxKernel(_type, var_name, params);
}

void AddMultiAuxVariableAction::addAuxVariable(std::string var_name)
{
	std::set<SubdomainID> blocks = getSubdomainIDs();

	if (_scalar_var)
		_problem->addAuxScalarVariable(var_name, _fe_type.order);

	else
	{
		if (_fe_type.order > 9)
			mooseError("Non-scalar AuxVariables must be CONSTANT, FIRST, SECOND, THIRD, FOURTH, FIFTH, SIXTH, SEVENTH, EIGHTH or NINTH order (" << _fe_type.order << " supplied)");

		if (blocks.empty())
			_problem->addAuxVariable(var_name, _fe_type);
		else
			_problem->addAuxVariable(var_name, _fe_type, &blocks);
	}
}
