
#include "CommonPostProcessorAction.h"
#include "MooseApp.h"
#include "FEProblem.h"
#include "MooseObjectAction.h"
#include "ActionFactory.h"
#include "Conversion.h"

template<>
InputParameters validParams<CommonPostProcessorAction>()
{
   InputParameters params = validParams<Action>();

   params.addParam<bool>("time_step", true, "Output the results using the default settings for Nemesis output");
   params.addParam<bool>("active_time", true, "Output the results using the default settings for Console output");
   params.addParam<bool>("residual", true, "Output the results using the default settings for XDA/XDR output (binary)");
   params.addParam<bool>("physic_time", false, " ");
   params.addParam<bool>("alive_time", false, "Output the scalar variable and postprocessors to a *.csv file using the default CSV output.");
   params.addParam<bool>("linear_step", false, "Output the results using the default settings for VTKOutput output");
   params.addParam<bool>("nonlinear_step", false, "Output the results using the default settings for XDA/XDR output (ascii)");

  return params;
}

CommonPostProcessorAction::CommonPostProcessorAction(const std::string & name, InputParameters params) :
    Action(name, params)
{
}

void CommonPostProcessorAction::act()
{
	if(getParam<bool>("alive_time"))
	{
		InputParameters params = _factory.getValidParams("RunTime");
		params.set<MooseEnum>("time_type") = "alive";
		_problem->addPostprocessor("RunTime", "alive_time", params);
	}

	if(getParam<bool>("active_time"))
	{
		InputParameters params = _factory.getValidParams("RunTime");
		params.set<MooseEnum>("time_type") = "active";
		_problem->addPostprocessor("RunTime", "active_time", params);
	}

	if(getParam<bool>("time_step"))
	{
		InputParameters params = _factory.getValidParams("NumTimeStep");
		_problem->addPostprocessor("NumTimeStep", "time_step", params);
	}

	if(getParam<bool>("nonlinear_step"))
	{
		InputParameters params = _factory.getValidParams("NumLinearIterations");
		_problem->addPostprocessor("NumLinearIterations", "linear_step", params);
	}

	if(getParam<bool>("linear_step"))
	{
		InputParameters params = _factory.getValidParams("NumNonlinearIterations");
		_problem->addPostprocessor("NumNonlinearIterations", "nonlinear_step", params);
	}

	if(getParam<bool>("residual"))
	{
		InputParameters params = _factory.getValidParams("VariableResidual");
//		params.set<VariableName>("variable") = "rho";
		_problem->addPostprocessor("VariableResidual", "residual", params);
	}


}

void CommonPostProcessorAction::create(std::string type)
{

}
