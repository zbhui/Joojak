#include "JoojakApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"

/// 单元积分
#include "EulerCellKernel.h"


/// 面积分
#include "EulerFaceKernel.h"

/// 初始条件
#include "IsoVortexIC.h"


/// 边界条件


/// 函数


/// 辅助kernel


/// Action


/// 材料属性
#include "EulerMaterial.h"

/// 时间步长增加策略

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

	/// 注册Kernel
	registerKernel(EulerCellKernel);

	/// 注册DGKernel
	registerDGKernel(EulerFaceKernel);

	/// 注册材料属性
	registerMaterial(EulerMaterial);
}

void
JoojakApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
}
