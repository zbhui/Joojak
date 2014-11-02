
#include "NumTimeStep.h"
#include "FEProblem.h"

template<>
InputParameters validParams<NumTimeStep>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  return params;
}

NumTimeStep::NumTimeStep(const std::string & name, InputParameters parameters) :
    GeneralPostprocessor(name, parameters),
    _feproblem(dynamic_cast<FEProblem &>(_subproblem))
{}

Real NumTimeStep::getValue()
{
  return _feproblem.timeStep();
}
