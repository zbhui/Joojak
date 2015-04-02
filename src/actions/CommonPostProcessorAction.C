
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

   params.addParam<bool>("physic_time", false, " ");
   params.addParam<bool>("time_step", false, "Output the results using the default settings for Nemesis output");
   params.addParam<bool>("active_time", true, "Output the results using the default settings for Console output");
   params.addParam<bool>("alive_time", false, "Output the scalar variable and postprocessors to a *.csv file using the default CSV output.");
   params.addParam<bool>("liner_ste", false, "Output the results using the default settings for VTKOutput output");
   params.addParam<bool>("nonlinear_step", false, "Output the results using the default settings for XDA/XDR output (ascii)");
   params.addParam<bool>("residual", false, "Output the results using the default settings for XDA/XDR output (binary)");
   params.addParam<unsigned int>("interval", 1, "The interval at which timesteps are output to the solution file");
   params.addParam<MultiMooseEnum>("output_on", Output::getExecuteOptions("timestep_end"), "Set to (initial|linear|nonlinear|timestep_end|timestep_begin|final|failed|custom) to execute only at that moment (default: timestep_end)");

  return params;
}

CommonPostProcessorAction::CommonPostProcessorAction(const std::string & name, InputParameters params) :
    Action(name, params),
	_action_params(_action_factory.getValidParams("AddPostprocessorAction"))
{
	  // Set the ActionWarehouse pointer in the parameters that will be passed to the actions created with this action
	  _action_params.set<ActionWarehouse *>("awh") = &_awh;
}

void CommonPostProcessorAction::act()
{
	if(getParam<bool>("alive_time"))
	{
		InputParameters params = _factory.getValidParams("RunTime");
		params.set<MooseEnum>("time_type") = "alive";
		_problem->addPostprocessor("RunTime", "run_time_new", params);

		params.set<MooseEnum>("time_type") = "active";
		_problem->addPostprocessor("RunTime", "active_time_new", params);
	}



}


