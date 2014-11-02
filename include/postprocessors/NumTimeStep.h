
#pragma once

#include "GeneralPostprocessor.h"

class NumTimeStep;

template<>
InputParameters validParams<NumTimeStep>();

class NumTimeStep: public GeneralPostprocessor
{
public:
	NumTimeStep(const std::string & name, InputParameters parameters);

  virtual void initialize() {}
  virtual void execute() {}

  virtual Real getValue();

protected:
  FEProblem & _feproblem;
};

