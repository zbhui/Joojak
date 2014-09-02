#include "JoojakApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

/// 单元积分
#include "EulerCellKernel.h"
#include "NSCellKernel.h"


/// 面积分
#include "EulerFaceKernel.h"
#include "NSFaceKernel.h"

/// 初始条件
#include "IsoVortexIC.h"
#include "CFDPassFlowIC.h"

/// 边界条件
#include "EulerBC.h"
#include "NSBC.h"

/// 函数
#include "IsoVortexExact.h"
#include "CouetteFlowExact.h"

/// 辅助kernel
#include "NSAuxVariable.h"

/// Action


/// 材料属性
#include "EulerCellMaterial.h"
#include "EulerFaceMaterial.h"
#include "EulerBndMaterial.h"
#include "IsoVortexBndMaterial.h"

#include "NSCellMaterial.h"
#include "NSFaceMaterial.h"
#include "NSBndMaterial.h"
#include "CouetteFlowBndMaterial.h"

/// 时间步长增加策略
#include "RatioTimeStepper.h"

/// PostProcessor
#include "CFDResidual.h"
#include "ElementExtremeTimeDerivative.h"
#include "CFDForcePostprocessor.h"

template<>
InputParameters validParams<JoojakApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

JoojakApp::JoojakApp(const std::string & name, InputParameters parameters) :
    MooseApp(name, parameters)
{
  srand(processor_id());

  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  JoojakApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  JoojakApp::associateSyntax(_syntax, _action_factory);
}

JoojakApp::~JoojakApp()
{
}

void
JoojakApp::registerApps()
{
  registerApp(JoojakApp);
}

void
JoojakApp::registerObjects(Factory & factory)
{
	/// 注册初始条件
	registerInitialCondition(IsoVortexIC);
	registerInitialCondition(CFDPassFlowIC);

	/// 注册边界条件
	registerBoundaryCondition(EulerBC);
	registerBoundaryCondition(NSBC);

	/// 注册Kernel
	registerKernel(EulerCellKernel);
	registerKernel(NSCellKernel);

	/// 注册DGKernel
	registerDGKernel(EulerFaceKernel);
	registerDGKernel(NSFaceKernel);

	/// 注册材料属性
	registerMaterial(EulerCellMaterial);
	registerMaterial(EulerFaceMaterial);
	registerMaterial(EulerBndMaterial);
	registerMaterial(IsoVortexBndMaterial);

	registerMaterial(NSCellMaterial);
	registerMaterial(NSFaceMaterial);
	registerMaterial(NSBndMaterial);
	registerMaterial(CouetteFlowBndMaterial);
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
}

void
JoojakApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
