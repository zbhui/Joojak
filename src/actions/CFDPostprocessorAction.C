
#include "CFDPostprocessorAction.h"
#include "MooseApp.h"
#include "FEProblem.h"

template<>
InputParameters validParams<CFDPostprocessorAction>()
{
  InputParameters params = validParams<Action>();
  MooseEnum time_options("alive active");
  params.addParam<MooseEnum>("time_type", time_options, "Whether to output the total elapsed or just the active time");
  return params;
}

CFDPostprocessorAction::CFDPostprocessorAction(const std::string & name, InputParameters params) :
    Action(name, params)
{
}

void CFDPostprocessorAction::act()
{
	InputParameters params = _factory.getValidParams("NumTimeStep");
	_problem->addPostprocessor("NumTimeStep", "num_timestep", params);

	params = _factory.getValidParams("Residual");
	_problem->addPostprocessor("Residual", "residual_final", params);

	params = _factory.getValidParams("CFDResidual");
	_problem->addPostprocessor("CFDResidual", "residual_initial", params);

	if(isParamValid("time_type"))
	{
		params = _factory.getValidParams("RunTime");
		params.set<MooseEnum>("time_type") = getParam<MooseEnum>("time_type");
		_problem->addPostprocessor("RunTime", "run_time", params);
	}

}

