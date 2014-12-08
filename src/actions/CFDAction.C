
#include "CFDAction.h"
#include "Factory.h"
#include "ActionFactory.h"
#include "MooseApp.h"
#include "MooseObject.h"
#include "Parser.h"
#include "FEProblem.h"

template<>
InputParameters validParams<CFDAction>()
{
  InputParameters params = validParams<Action>();
  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "The names of the convection and diffusion variables in the simulation");

  return params;
}

CFDAction::CFDAction(const std::string & name, InputParameters params) :
    Action(name, params)
{
}

void
CFDAction::act()
{
  std::vector<NonlinearVariableName> variables = getParam<std::vector<NonlinearVariableName> > ("variables");
  std::vector<VariableName> vel_vec_variable;

  /**
   * We need to manually setup our Convection-Diffusion and Diffusion variables on our two
   * variables we are expecting from the input file.  Much of the syntax below is hidden by the
   * parser system but we have to set things up ourselves this time.
   */

  // Do some error checking
  mooseAssert(variables.size() == 6, "Expected 2 variables, received " + variables.size());

  // Setup our Diffusion Kernel on the "u" variable
  {

//	InputParameters  _moose_object_par = validParams<MooseObject>();
//	_moose_object_par.get<Real>
//	MooseObject moose_object(name, parameters);
//	MooseObject moose_object(name(), parameters());
//	Real mach = moose_object.getParam<Real>("mach");
//	Real mach = _app.getParam<Real>("mach");
//	InputParameters pp = _action_factory.getValidParams("add_mesh");
    InputParameters params = _factory.getValidParams("CFDPressure");
    params.set<AuxVariableName>("variable") = variables[0];
    params.set<unsigned int>("component") = -1;
//    params.set<VariableName>("conservation_vale") = ""
    _problem->addAuxKernel("CFDPressure", "aux_pressure", params);
  }

/*  // Setup our Convection Kernel on the "u" variable coupled to the diffusion variable "v"
  {
    InputParameters params = _factory.getValidParams("Convection");
    params.set<NonlinearVariableName>("variable") = variables[0];
//    params.addCoupledVar("some_variable", "The gradient of this var");
    vel_vec_variable.push_back(variables[1]);
    params.set<std::vector<VariableName> >("some_variable") = vel_vec_variable;
    _problem->addKernel("Convection", "conv", params);
  }

  // Setup out Diffusion Kernel on the "v" variable
  {
    InputParameters params = _factory.getValidParams("Diffusion");
    params.set<NonlinearVariableName>("variable") = variables[1];
    _problem->addKernel("Diffusion", "diff_v", params);
  }*/
}
