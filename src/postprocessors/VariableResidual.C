
#include "VariableResidual.h"

#include "FEProblem.h"
#include "SubProblem.h"

template<>
InputParameters validParams<VariableResidual>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  return params;
}

VariableResidual::VariableResidual(const std::string & name, InputParameters parameters) :
    GeneralPostprocessor(name, parameters)
{}

Real VariableResidual::getValue()
{
  return _fe_problem.getNonlinearSystem()._initial_residual_before_preset_bcs;
//	TransientNonlinearImplicitSystem & _sys = _fe_problem.getNonlinearSystem().sys();
//	return _sys.calculate_norm(*_sys.rhs, 0, DISCRETE_L2);
}

