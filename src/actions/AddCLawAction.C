

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
  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "非线性变量");
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
		addVariable();
	else if(_current_task == "add_kernel")
		addKernel();
	else if(_current_task == "add_dg_kernel")
		addDGKernel();
	else if(_current_task == "add_aux_variable")
		addAuxVariable();
	else if(_current_task == "add_aux_kernel")
		addAuxKernel();
	else if(_current_task == "add_bc")
		addBoundaryCondition();
	else
	{
		mooseError("未定义task: " << _current_task);
	}
}

void AddCLawAction::addVariable()
{
	Real scale_factor = isParamValid("scaling") ? getParam<Real>("scaling") : 1;
	FEType fe_type(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
             Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family")));
	for (int i = 0; i < _variables.size(); ++i)
	{
		_problem->addVariable(_variables[i], fe_type, scale_factor);
	}
}

void AddCLawAction::addKernel()
{
	std::string time_kernel_name = "TimeDerivative";
	InputParameters params = _factory.getValidParams(time_kernel_name);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addKernel(time_kernel_name, _variables[i] + "_time", params);
	}

	std::string cell_kernel_name = "CLawCellKernel";
	params = _factory.getValidParams(cell_kernel_name);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addKernel(cell_kernel_name, _variables[i] + "_space", params);
	}
}

void AddCLawAction::addDGKernel()
{
	std::string face_kernel_name = "CLawFaceKernel";
	InputParameters params = _factory.getValidParams(face_kernel_name);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addDGKernel(face_kernel_name, _variables[i] + "_dg", params);
	}
}

void AddCLawAction::addAuxVariable()
{
	FEType fe_type(Utility::string_to_enum<Order>("CONSTANT"), Utility::string_to_enum<FEFamily>("MONOMIAL"));

	_problem->addAuxVariable("proc_id", fe_type);
}

void AddCLawAction::addAuxKernel()
{
	InputParameters params = _factory.getValidParams("ProcessorIDAux");
	params.set<AuxVariableName>("variable") = "proc_id";
	_problem->addAuxKernel("ProcessorIDAux", "proc_id", params);
}

void AddCLawAction::addBoundaryCondition()
{
    std::vector<BoundaryName> boundary;
    std::set<boundary_id_type> boundary_id =  _mesh->meshBoundaryIds();
    for (std::set<boundary_id_type>::const_iterator itor = boundary_id.begin(); itor != boundary_id.end(); ++itor)
    	boundary.push_back(_mesh->getMesh().get_boundary_info().sideset_name(*itor));

	std::string boun_cond_name = "NSBC";
    InputParameters params = _factory.getValidParams(boun_cond_name);
    params.set<std::vector<BoundaryName> >("boundary") = boundary;

	for (int i = 0; i < _variables.size(); ++i)
	{
	    params.set<NonlinearVariableName>("variable") = _variables[i];
	    _problem->addBoundaryCondition(boun_cond_name, _variables[i]+"_bc", params);
	}
}

