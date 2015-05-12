
#include "Moose.h"
#include "JoojakApp.h"

#include "CLawICAction.h"
#include "CLawAuxVariablesAction.h"
#include "JoojakRevision.h"
#include "AppFactory.h"
#include "ActionFactory.h"
#include "Syntax.h"

/// Action
#include "CFDAddVariablesAction.h"
#include "CFDBoundaryConditionAction.h"
#include "CFDKernelsAction.h"
#include "CFDDGKernelsAction.h"
#include "CFDPostprocessorAction.h"
#include "CommonPostProcessorAction.h"
#include "AddCLawAction.h"
#include "AddMultiVariableAction.h"
#include "AddMultiAuxVariableAction.h"

/// 单元积分
#include "CLawCellKernel.h"
#include "EulerCellKernel.h"
#include "NSCellKernel.h"
#include "KOCellKernel.h"

#include "EmptyTimeDerivative.h"
#include "ElasticityKernel.h"

/// 面积分
#include "CLawFaceKernel.h"
#include "EulerFaceKernel.h"
#include "NSFaceKernel.h"
#include "KOFaceKernel.h"

/// 初始条件
#include "IsoVortexIC.h"
#include "CFDPassFlowIC.h"
#include "KOIC.h"

/// 边界条件
#include "EulerBC.h"
#include "NSBC.h"
#include "KOBC.h"
#include "CLawBoundaryCondition.h"

/// 函数
#include "IsoVortexExact.h"
#include "CouetteFlowExact.h"

/// 辅助kernel
#include "NSAuxVariable.h"
#include "EmptyTimeDerivative.h"

/// 材料属性
#include "EulerCellMaterial.h"
#include "EulerFaceMaterial.h"
#include "EulerBndMaterial.h"
#include "IsoVortexBndMaterial.h"

#include "NSCellMaterial.h"
#include "NSFaceMaterial.h"
#include "NSBndMaterial.h"
#include "CouetteFlowBndMaterial.h"

#include "CLawFaceMaterial.h"
#include "CLawCellMaterial.h"
#include "CLawBoundaryMaterial.h"
#include "KOCellMaterial.h"
#include "KOFaceMaterial.h"
#include "KOBndMaterial.h"

#include "LinearElasticityMaterial.h"
/// 时间步长增加策略
#include "RatioTimeStepper.h"

/// PostProcessor
#include "CFDResidual.h"
#include "ElementExtremeTimeDerivative.h"
#include "CFDForcePostprocessor.h"
#include "BumpElementL2Error.h"
#include "NumTimeStep.h"
#include "VariableResidual.h"
#include "IsoVortexElementL2Error.h"

/// UserObject
#include "CFDForceUserObject.h"

/// VectorPostProcessor
#include "PressureAndSkinFrictionCoeff.h"

/// Executioner
#include "SteadyTransientExecutioner.h"

#include "SAInclude.h"

/// mesh modifier
#include "BuildSideSetFromBlock.h"

/// problem
#include "CLawProblem.h"
#include "NavierStokesProblem.h"
#include "EulerProblem.h"
#include "SAProblem.h"

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
	registerInitialCondition(IsoVortexIC);
	registerInitialCondition(CFDPassFlowIC);
	registerInitialCondition(KOIC);

	/// 注册边界条件
	registerBoundaryCondition(EulerBC);
	registerBoundaryCondition(NSBC);
	registerBoundaryCondition(KOBC);
	registerBoundaryCondition(CLawBoundaryCondition);

	/// 注册Kernel
	registerKernel(CLawCellKernel);
	registerKernel(EulerCellKernel);
	registerKernel(NSCellKernel);
	registerKernel(KOCellKernel);

	registerKernel(EmptyTimeDerivative);
	registerKernel(ElasticityKernel);

	/// 注册DGKernel
	registerDGKernel(CLawFaceKernel);
	registerDGKernel(EulerFaceKernel);
	registerDGKernel(NSFaceKernel);
	registerDGKernel(KOFaceKernel);

	/// 注册材料属性
	registerMaterial(EulerCellMaterial);
	registerMaterial(EulerFaceMaterial);
	registerMaterial(EulerBndMaterial);
	registerMaterial(IsoVortexBndMaterial);

	registerMaterial(NSCellMaterial);
	registerMaterial(NSFaceMaterial);
	registerMaterial(NSBndMaterial);
	registerMaterial(CouetteFlowBndMaterial);

	registerMaterial(KOCellMaterial);
	registerMaterial(KOFaceMaterial);
	registerMaterial(KOBndMaterial);

    registerMaterial(LinearElasticityMaterial);

    registerMaterial(CLawFaceMaterial);
    registerMaterial(CLawCellMaterial);
    registerMaterial(CLawBoundaryMaterial);
	/// 注册函数
	registerFunction(IsoVortexExact);
	registerFunction(CouetteFlowExact);

	///注册辅助kernel
	registerAux(NSAuxVariable);

	/// 注册时间步长
	registerTimeStepper(RatioTimeStepper);

	/// 注册后处理
	registerPostprocessor(CFDResidual);
	registerPostprocessor(ElementExtremeTimeDerivative);
	registerPostprocessor(CFDForcePostprocessor);
	registerPostprocessor(BumpElementL2Error);
	registerPostprocessor(NumTimeStep);
	registerPostprocessor(VariableResidual);
	registerPostprocessor(IsoVortexElementL2Error);


	/// 注册UserObject
	registerPostprocessor(CFDForceUserObject);


	registerVectorPostprocessor(PressureAndSkinFrictionCoeff);

	registerExecutioner(SteadyTransientExecutioner);

	registerMeshModifier(BuildSideSetFromBlock);

	registerProblem(CLawProblem);
	registerProblem(NavierStokesProblem);
	registerProblem(EulerProblem);
	registerProblem(SAProblem);

	registerSAObjects(factory);
}

void JoojakApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
	/// 注册Action
	syntax.registerActionSyntax("CLawAuxVariablesAction", "AuxVariables");
	syntax.registerActionSyntax("CLawICAction", "ICs", "add_ic");
	syntax.registerActionSyntax("CommonPostProcessorAction", "Postprocessors", "add_postprocessor");
//	syntax.registerActionSyntax("AddCLawAction", "Problem");

	registerAction(CLawICAction, "add_ic");
	registerAction(CFDPostprocessorAction, "add_postprocessor");
	registerAction(CLawAuxVariablesAction, "add_aux_variable");
	registerAction(CLawAuxVariablesAction, "add_aux_kernel");

	registerAction(CommonPostProcessorAction, "add_postprocessor");
//	registerAction(AddCLawAction, "add_variable");
//	registerAction(AddCLawAction, "add_kernel");
//	registerAction(AddCLawAction, "add_dg_kernel");
//	registerAction(AddCLawAction, "add_aux_variable");
//	registerAction(AddCLawAction, "add_aux_kernel");
//	registerAction(AddCLawAction, "add_bc");

	syntax.registerActionSyntax("AddMultiVariableAction", "Problem/Variables");
	registerAction(AddMultiVariableAction, "add_variable");
	registerAction(AddMultiVariableAction, "add_kernel");
	registerAction(AddMultiVariableAction, "add_dg_kernel");
	registerAction(AddMultiVariableAction, "add_bc");

	syntax.registerActionSyntax("AddMultiAuxVariableAction", "Problem/AuxVariables/*");
	registerAction(AddMultiAuxVariableAction, "add_aux_variable");
	registerAction(AddMultiAuxVariableAction, "add_aux_kernel");



}

void JoojakApp::registerSAObjects(Factory & factory)
{
	/// 注册初始条件
	registerInitialCondition(SAIC);

	/// 注册边界条件
	registerBoundaryCondition(SABC);

	/// 注册Kernel
	registerKernel(SACellKernel);

	/// 注册DGKernel
	registerDGKernel(SAFaceKernel);

	/// 注册材料属性
	registerMaterial(SACellMaterial);
	registerMaterial(SAFaceMaterial);
	registerMaterial(SABndMaterial);

	/// 注册函数

	///注册辅助kernel
	registerAux(SAAuxVariable);
	registerAux(NearestWallDistance);

	/// 注册时间步长

	/// 注册后处理
}
