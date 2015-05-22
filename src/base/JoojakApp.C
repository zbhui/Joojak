
#include "Moose.h"
#include "JoojakApp.h"

#include "JoojakRevision.h"
#include "AppFactory.h"
#include "ActionFactory.h"
#include "Syntax.h"

/// Action
#include "CFDAuxVariable.h"
#include "CLawAuxVariablesAction.h"
#include "CLawICAction.h"
#include "CommonPostProcessorAction.h"
#include "AddMultiVariableAction.h"
#include "AddMultiAuxVariableAction.h"

/// 单元积分
#include "CLawCellKernel.h"

#include "EmptyTimeDerivative.h"
#include "ElasticityKernel.h"

/// 面积分
#include "CLawFaceKernel.h"

/// 初始条件
#include "CLawIC.h"
#include "CFDPassFlowIC.h"

/// 边界条件
#include "CLawBoundaryCondition.h"

/// 函数
#include "CouetteFlowExact.h"

/// 辅助kernel
#include "EmptyTimeDerivative.h"

/// 材料属性
#include "CLawFaceMaterial.h"
#include "CLawCellMaterial.h"
#include "CLawBoundaryMaterial.h"

#include "LinearElasticityMaterial.h"
/// 时间步长增加策略
#include "RatioTimeStepper.h"

/// PostProcessor
#include "CFDResidual.h"
#include "ElementExtremeTimeDerivative.h"
#include "NumTimeStep.h"
#include "VariableResidual.h"
#include "IsoVortexElementL2Error.h"
#include "CouetteFlowElementL2Error.h"

/// UserObject
#include "CFDForceUserObject.h"

/// VectorPostProcessor

/// Executioner
#include "SteadyTransientExecutioner.h"

/// mesh modifier
#include "BuildSideSetFromBlock.h"

/// problem
#include "CLawProblem.h"
#include "NavierStokesProblem.h"
#include "EulerProblem.h"
#include "SAProblem.h"
#include "IsoVortexProblem.h"
#include "CouetteFlowProblem.h"
#include "SodProblem.h"

template<>
InputParameters validParams<JoojakApp>()
{
  InputParameters params = validParams<MooseApp>();
  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

JoojakApp::JoojakApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  JoojakApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  JoojakApp::associateSyntax(_syntax, _action_factory);
}

JoojakApp::~JoojakApp()
{
}

void JoojakApp::printHeader()
{
	std::string line("*********************************\n\n");
	Moose::out << COLOR_CYAN << line << COLOR_DEFAULT;
	Moose::out << "计算流体力学间断有限元计算器 JOOJAK \n\n";
	Moose::out << "Joojak version: " <<  COLOR_MAGENTA << JOOJAK_REVISION << COLOR_DEFAULT << std::endl << std::endl;
	Moose::out << COLOR_CYAN << line << COLOR_DEFAULT;

//
//	std::string joojak =
//			std::string("		     _             _       _\n")+
//			std::string("		    | | ___   ___ (_) __ _| | __\n")+
//			std::string("		 _  | |/ _ \ / _ \| |/ _` | |/ /\n")+
//			std::string("		| |_| | (_) | (_) | | (_| |   <\n")+
//			std::string("		 \___/ \___/ \___// |\__,_|_|\_\\n")+
//			std::string("		                |__/\n");
//
//Moose::out << joojak << std::endl;
}

void JoojakApp::run()
{
	printHeader();
	setupOptions();
	runInputFile();
	executeExecutioner();
}
void JoojakApp::registerApps()
{
  registerApp(JoojakApp);
}

void JoojakApp::registerObjects(Factory & factory)
{
	/// 注册初始条件
	registerInitialCondition(CFDPassFlowIC);
	registerInitialCondition(CLawIC);

	/// 注册边界条件
	registerBoundaryCondition(CLawBoundaryCondition);

	/// 注册Kernel
	registerKernel(CLawCellKernel);

	registerKernel(EmptyTimeDerivative);
	registerKernel(ElasticityKernel);

	/// 注册DGKernel
	registerDGKernel(CLawFaceKernel);

	/// 注册材料属性
    registerMaterial(LinearElasticityMaterial);

    registerMaterial(CLawFaceMaterial);
    registerMaterial(CLawCellMaterial);
    registerMaterial(CLawBoundaryMaterial);
	/// 注册函数
	registerFunction(CouetteFlowExact);

	///注册辅助kernel
	registerAux(CFDAuxVariable);

	/// 注册时间步长
	registerTimeStepper(RatioTimeStepper);

	/// 注册后处理
	registerPostprocessor(CFDResidual);
	registerPostprocessor(ElementExtremeTimeDerivative);
	registerPostprocessor(NumTimeStep);
	registerPostprocessor(VariableResidual);
	registerPostprocessor(IsoVortexElementL2Error);
	registerPostprocessor(CouetteFlowElementL2Error);


	/// 注册UserObject
	registerPostprocessor(CFDForceUserObject);

	registerExecutioner(SteadyTransientExecutioner);

	registerMeshModifier(BuildSideSetFromBlock);


	/// 注册Problem
//	registerProblem(CLawProblem);
	registerProblem(NavierStokesProblem);
	registerProblem(EulerProblem);
	registerProblem(SAProblem);
	registerProblem(IsoVortexProblem);
	registerProblem(CouetteFlowProblem);
	registerProblem(SodProblem);
}

void JoojakApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
	/// 注册Action

	syntax.registerActionSyntax("CLawICAction", "ICs", "add_ic");
	registerAction(CLawICAction, "add_ic");

	syntax.registerActionSyntax("CLawAuxVariablesAction", "AuxVariables");
	registerAction(CLawAuxVariablesAction, "add_aux_variable");
	registerAction(CLawAuxVariablesAction, "add_aux_kernel");

	syntax.registerActionSyntax("CommonPostProcessorAction", "Postprocessors", "add_postprocessor");
	registerAction(CommonPostProcessorAction, "add_postprocessor");

	syntax.registerActionSyntax("AddMultiVariableAction", "Problem/Variables");
	registerAction(AddMultiVariableAction, "add_variable");
	registerAction(AddMultiVariableAction, "add_kernel");
	registerAction(AddMultiVariableAction, "add_dg_kernel");
	registerAction(AddMultiVariableAction, "add_bc");

	syntax.registerActionSyntax("AddMultiAuxVariableAction", "Problem/AuxVariables/*");
	registerAction(AddMultiAuxVariableAction, "add_aux_variable");
	registerAction(AddMultiAuxVariableAction, "add_aux_kernel");

}

