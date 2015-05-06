
#include <sstream>
#include <stdexcept>

#include "AddMultiVariableAction.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "EigenSystem.h"

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"

template<>
InputParameters validParams<AddMultiVariableAction>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());
  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "非线性变量");

  return params;
}

AddMultiVariableAction::AddMultiVariableAction(const std::string & name, InputParameters params) :
	AddVariableAction(name, params),
	_variables(getParam<std::vector<NonlinearVariableName> >("variables"))

{
}

void AddMultiVariableAction::act()
{
	for (int i = 0; i < _variables.size(); ++i)
	{
		std::string var_name = _variables[i];
		if (_current_task == "add_variable")
			addVariable(var_name);
		if (_current_task == "add_kernel")
			addKernel(var_name, i);
		if (_current_task == "add_dg_kernel")
			addDGKernel(var_name, i);
		if (_current_task == "add_bc")
			addBoundaryCondition(var_name, i);
		if (_current_task == "add_ic")
			setInitialCondition(var_name, i);
	}
}

void AddMultiVariableAction::addKernel(std::string var_name, int i)
{
	std::string time_kernel_name = "TimeDerivative";
	InputParameters params = _factory.getValidParams(time_kernel_name);
	params.set<NonlinearVariableName>("variable") = var_name;
	_problem->addKernel(time_kernel_name, var_name + "_time", params);

	std::string cell_kernel_name = "CLawCellKernel";
	params = _factory.getValidParams(cell_kernel_name);
	params.set<NonlinearVariableName>("variable") = var_name;
	params.set<int>("component") = i;
	_problem->addKernel(cell_kernel_name, var_name + "_space", params);
}

void AddMultiVariableAction::addDGKernel(std::string var_name, int i)
{
	std::string face_kernel_name = "CLawFaceKernel";
	InputParameters params = _factory.getValidParams(face_kernel_name);
	params.set<NonlinearVariableName>("variable") = _variables[i];
	params.set<int>("component") = i;
	_problem->addDGKernel(face_kernel_name,  var_name + "_dg", params);
}

void AddMultiVariableAction::addBoundaryCondition(std::string var_name, int i)
{
    std::vector<BoundaryName> boundary;
    std::set<boundary_id_type> boundary_id =  _mesh->meshBoundaryIds();
    for (std::set<boundary_id_type>::const_iterator itor = boundary_id.begin(); itor != boundary_id.end(); ++itor)
    	boundary.push_back(_mesh->getMesh().get_boundary_info().sideset_name(*itor));

    std::string boun_cond_name = "CLawBoundaryCondition";
    InputParameters params = _factory.getValidParams(boun_cond_name);
    params.set<std::vector<BoundaryName> >("boundary") = boundary;

    params.set<NonlinearVariableName>("variable") = var_name;
    params.set<int>("component") = i;
    _problem->addBoundaryCondition(boun_cond_name, var_name +"_bc", params);
}

void AddMultiVariableAction::setInitialCondition(std::string var_name, int eq)
{
	Real initial = getParam<Real>("initial_condition");
	if (initial > _abs_zero_tol || initial < -_abs_zero_tol)
	{
		if (_scalar_var)
		{
			// built a ScalarConstantIC object
			InputParameters params = _factory.getValidParams("ScalarConstantIC");
			params.set<VariableName>("variable") = var_name;
			params.set<Real>("value") = initial;
			_problem->addInitialCondition("ScalarConstantIC", "ic", params);
		}
		else
		{
			// built a ConstantIC object
			InputParameters params = _factory.getValidParams("ConstantIC");
			params.set<VariableName>("variable") = var_name;
			params.set<Real>("value") = initial;
			_problem->addInitialCondition("ConstantIC", "ic", params);
		}
	}
}

