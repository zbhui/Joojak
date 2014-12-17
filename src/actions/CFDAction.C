
#include "CFDAction.h"
#include "AddVariableAction.h"
#include "Factory.h"
#include "ActionFactory.h"
#include "MooseApp.h"
#include "MooseObject.h"
#include "Parser.h"
#include "FEProblem.h"

#include "libmesh/string_to_enum.h"

template<>
InputParameters validParams<CFDAction>()
{
  MooseEnum families(AddVariableAction::getNonlinearVariableFamilies());
  MooseEnum orders(AddVariableAction::getNonlinearVariableOrders());
  MooseEnum equations("EULER NS SA", "EULER");
  MooseEnum init_cond("CFDPassFlowIC IsoVortexIC", "CFDPassFlowIC");
  InputParameters params = validParams<Action>();
  params.addParam<MooseEnum>("family", families, "Specifies the family of FE shape functions to use for this variable");
  params.addParam<MooseEnum>("order", orders,  "Specifies the order of the FE shape function to use for this variable (additional orders not listed are allowed)");
  params.addRequiredParam<MooseEnum>("type", equations, "求解方程组类型");
  params.addRequiredParam<MooseEnum>("init_cond", init_cond, "求解方程组类型");
  params.addParam<std::vector<BoundaryName> >("boundary", "The list of boundary IDs from the mesh where the object will be applied");
  params.addParam<std::vector<AuxVariableName> >("aux_variables", "辅助变量名");
  return params;
}

CFDAction::CFDAction(const std::string & name, InputParameters params) :
    Action(name, params),
	_aux_variables(getParam<std::vector<AuxVariableName> >("aux_variables"))
{
	MooseEnum type(getParam<MooseEnum>("type"));
	MooseEnum init_cond(getParam<MooseEnum>("init_cond"));
	if(type == "EULER")
	{
		_variables.push_back("rho");
		_variables.push_back("momentum_x");
		_variables.push_back("momentum_y");
		_variables.push_back("momentum_z");
		_variables.push_back("rhoe");

		_aux_kernel_name = "NSAuxVariable";
		_time_kernel_name = "TimeDerivative";
		_cell_kernel_name = "EulerCellKernel";
		_face_kernel_name = "EulerFaceKernel";
		_cell_material_name = "EulerCellMaterial";
		_face_material_name = "EulerFaceMaterial";
		_init_cond_name = std::string(init_cond);
		_boun_cond_name = "EulerBC";

		if(_init_cond_name != "CFDPassFlowIC" && _init_cond_name != "IsoVortexIC")
			mooseError("Euler方程没有这种初始条件 " << _init_cond_name);
	}
	else if(type == "NS")
	{
		_variables.push_back("rho");
		_variables.push_back("momentum_x");
		_variables.push_back("momentum_y");
		_variables.push_back("momentum_z");
		_variables.push_back("rhoe");
	}
	else if(type == "SA")
	{
		_variables.push_back("rho");
		_variables.push_back("momentum_x");
		_variables.push_back("momentum_y");
		_variables.push_back("momentum_z");
		_variables.push_back("rhoe");
		_variables.push_back("rhon");
	}
	else
	{
		mooseError("未知的方程类型");
	}


}

void CFDAction::act()
{
	if (_current_task == "add_variable")
	{
		addVariables();
	}
	else if (_current_task == "add_aux_variable")
	{
		addAuxVariables();
	}
	else if (_current_task == "add_ic")
	{
		setInitialCondition();
	}
	else if (_current_task == "add_bc")
	{
		setBoundaryCondition();
	}
	else if (_current_task == "add_kernel")
	{
		addKernel();
	}
	else if (_current_task == "add_aux_kernel")
	{
		addAuxKernel();
	}
	else if (_current_task == "add_dg_kernel")
	{
		addDGKernel();
	}
	else if (_current_task == "add_material")
	{
//		addMaterial();
	}
	else
	{

	}

}

void CFDAction::addVariables()
{
	Real scale_factor = isParamValid("scaling") ? getParam<Real>("scaling") : 1;
	FEType fe_type(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
             Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family")));
	for (int i = 0; i < _variables.size(); ++i)
	{
		_problem->addVariable(_variables[i], fe_type, scale_factor);
	}
}

void CFDAction::addAuxVariables()
{
	FEType fe_type(Utility::string_to_enum<Order>(getParam<MooseEnum>("order")),
             Utility::string_to_enum<FEFamily>(getParam<MooseEnum>("family")));

	for (int i = 0; i < _aux_variables.size(); ++i)
	{
		_problem->addAuxVariable(_aux_variables[i], fe_type);
	}
}

void CFDAction::setInitialCondition()
{
    InputParameters params = _factory.getValidParams(_init_cond_name);
	for (int i = 0; i < _variables.size(); ++i)
	{
	    params.set<VariableName>("variable") = _variables[i];
	    _problem->addInitialCondition(_init_cond_name, _variables[i]+"_ic", params);
	}
}

void CFDAction::setBoundaryCondition()
{
    InputParameters params = _factory.getValidParams("EulerBC");
    params.set<std::vector<BoundaryName> >("boundary") = getParam<std::vector<BoundaryName> >("boundary");
	for (int i = 0; i < _variables.size(); ++i)
	{
	    params.set<NonlinearVariableName>("variable") = _variables[i];
	    _problem->addBoundaryCondition("EulerBC", _variables[i]+"_bc", params);
	}
}

void CFDAction::addKernel()
{
	InputParameters params = _factory.getValidParams(_time_kernel_name);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addKernel(_time_kernel_name, _variables[i] + "_time", params);
	}

	params = _factory.getValidParams(_cell_kernel_name);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addKernel(_cell_kernel_name, _variables[i] + "_space", params);
	}
}

void CFDAction::addAuxKernel()
{
	std::vector<VariableName> var_name;
	for (int i = 0; i < _variables.size(); ++i)
		var_name.push_back(_variables[i]);

	InputParameters params = _factory.getValidParams(_aux_kernel_name);
	params.set<std::vector<VariableName> >("variables") = var_name;

	for (int i = 0; i < _aux_variables.size(); ++i)
	{
		params.set<AuxVariableName>("variable") = _aux_variables[i];
		_problem->addAuxKernel(_aux_kernel_name, _aux_variables[i], params);
	}
}

void CFDAction::addDGKernel()
{
    InputParameters params = _factory.getValidParams(_face_kernel_name);
	for (int i = 0; i < _variables.size(); ++i)
	{
		params.set<NonlinearVariableName>("variable") = _variables[i];
		_problem->addDGKernel(_face_kernel_name, _variables[i] + "_dg", params);
	}
}
void CFDAction::addMaterial()
{
	std::vector<VariableName> var_name;
	for (int i = 0; i < _variables.size(); ++i)
		var_name.push_back(_variables[i]);

	InputParameters params = _factory.getValidParams("EulerCellMaterial");
	params.set<std::vector<VariableName> >("variables") = var_name;
	params.set<std::vector<SubdomainName> >("block") = getParam<std::vector<SubdomainName> >("block");
	_problem->addMaterial("EulerCellMaterial", "Materials/cell_mateiral1" , params);
	_problem->addMaterial("EulerFaceMaterial", "Materials/face_mateiral1" , params);

	params = _factory.getValidParams("IsoVortexBndMaterial");
//	_app.parser().extractParams(_name, params);
	params.set<std::vector<VariableName> >("variables") = var_name;
	params.set<std::vector<BoundaryName> >("boundary") = getParam<std::vector<BoundaryName> >("boundary");
	_problem->addMaterial("IsoVortexBndMaterial", "bnd_mateiral" , params);


}
