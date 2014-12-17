#include "CFDMetaAction.h"
#include "Factory.h"
#include "ActionFactory.h"
#include "MooseObjectAction.h"
#include "MooseApp.h"
#include "FEProblem.h"

#include "libmesh/vector_value.h"

template<>
InputParameters validParams<CFDMetaAction>()
{
	InputParameters params = validParams<Action>();
	params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "CFD variables name");
    params.addRequiredParam<std::vector<BoundaryName> >("boundary", "The list of boundary IDs from the mesh where the pressure will be applied");

	return params;
}

CFDMetaAction::CFDMetaAction(const std::string & name, InputParameters params) :
		Action(name, params),
		_boundary(getParam<std::vector<BoundaryName> >("boundary"))
{
}

void CFDMetaAction::act()
{
	MooseSharedPointer<Action> action;
	MooseSharedPointer<MooseObjectAction> moose_object_action;

	std::vector<NonlinearVariableName> variables = getParam<std::vector<NonlinearVariableName> > ("variables");

	/**
	 * We need to manually setup our Convection-Diffusion and Diffusion variables on our two
	 * variables we are expecting from the input file.  Much of the syntax below is hidden by the
	 * parser system but we have to set things up ourselves this time.
	 */

	// Do some error checking
	mooseAssert(variables.size() == 5, "Expected 5 variables");

	//*******************************************//
	//**************** Variables ****************//
	//*******************************************//
	InputParameters action_params = _action_factory.getValidParams("AddVariableAction");
//	action_params.set<ActionWarehouse *>("awh") = &_awh;
//	for (unsigned int i=0; i<variables.size(); ++i)
//	{
//		action = _action_factory.create("AddVariableAction", "Variables/" + variables[i], action_params);
//		_awh.addActionBlock(action);
//	}

	//*******************************************//
	//****************    ICs    ****************//
	//*******************************************//
//	action_params = _action_factory.getValidParams("AddICAction");
//	action_params.set<ActionWarehouse *>("awh") = &_awh;
//	action_params.set<std::string>("type") = "IsoVortexIC";
//	for (unsigned int i=0; i<variables.size(); ++i)
//	{
//		action = _action_factory.create("AddICAction", "ICs/" + variables[i]+"_ic", action_params);
//		moose_object_action = MooseSharedNamespace::dynamic_pointer_cast<MooseObjectAction>(action);
//		mooseAssert (moose_object_action.get(), "Dynamic Cast failed");
//		InputParameters & params = moose_object_action->getObjectParams();
//		params.set<VariableName>("variable") = variables[i];
//		_awh.addActionBlock(action);
//	}


	//*******************************************//
	//**************** Kernels ******************//
	//*******************************************//
	action_params = _action_factory.getValidParams("AddKernelAction");
	action_params.set<ActionWarehouse *>("awh") = &_awh;
	action_params.set<std::string>("type") = "EulerCellKernel";
	for (unsigned int i=0; i<variables.size(); ++i)
	{
		action = _action_factory.create("AddKernelAction", "Kernels/"+variables[i]+"_space", action_params);
		moose_object_action = MooseSharedNamespace::dynamic_pointer_cast<MooseObjectAction>(action);
		mooseAssert (moose_object_action.get(), "Dynamic Cast failed");

		InputParameters & params = moose_object_action->getObjectParams();
		params.set<NonlinearVariableName>("variable") = variables[i];
		_awh.addActionBlock(action);

	}

	action_params = _action_factory.getValidParams("AddKernelAction");
	action_params.set<ActionWarehouse *>("awh") = &_awh;
	action_params.set<std::string>("type") = "TimeDerivative";
	for (unsigned int i=0; i<variables.size(); ++i)
	{
		action = _action_factory.create("AddKernelAction", "Kernels/"+variables[i]+"_time", action_params);
		moose_object_action = MooseSharedNamespace::dynamic_pointer_cast<MooseObjectAction>(action);
		mooseAssert (moose_object_action.get(), "Dynamic Cast failed");

		InputParameters & params = moose_object_action->getObjectParams();
		params.set<NonlinearVariableName>("variable") = variables[i];
		_awh.addActionBlock(action);
	}

	//*******************************************//
	//**************** DGKernels ****************//
	//*******************************************//
	action_params = _action_factory.getValidParams("AddDGKernelAction");
	action_params.set<ActionWarehouse *>("awh") = &_awh;
	action_params.set<std::string>("type") = "EulerFaceKernel";
	for (unsigned int i=0; i<variables.size(); ++i)
	{
		action = _action_factory.create("AddDGKernelAction", "DGKernels/"+variables[i]+"_dg", action_params);
		moose_object_action = MooseSharedNamespace::dynamic_pointer_cast<MooseObjectAction>(action);
		mooseAssert (moose_object_action.get(), "Dynamic Cast failed");

		InputParameters & params = moose_object_action->getObjectParams();
		params.set<NonlinearVariableName>("variable") = variables[i];
		_awh.addActionBlock(action);
	}


	//*******************************************//
	//****************    BCs    ****************//
	//*******************************************//
////	int i = 0;
	action_params = _action_factory.getValidParams("AddBCAction");
	action_params.set<ActionWarehouse *>("awh") = &_awh;

		action_params.set<std::string>("type") = "EulerBC";
	for (unsigned int i=0; i<variables.size(); ++i)
	{
		action_params.set<NonlinearVariableName>("variable") = variables[i];
		action_params.set<std::vector<BoundaryName> >("boundary") = _boundary;
		action = _action_factory.create("AddBCAction", "BCs/"+variables[i]+"_bc", action_params);
		moose_object_action = MooseSharedNamespace::dynamic_pointer_cast<MooseObjectAction>(action);
		mooseAssert (moose_object_action.get(), "Dynamic Cast failed");

		InputParameters & params = moose_object_action->getObjectParams();
		params.set<std::vector<BoundaryName> >("boundary") = _boundary;
		params.set<NonlinearVariableName>("variable") = variables[i];
		params += validParams<GlobalParamsAction>();
		_awh.addActionBlock(action);
	}


}
