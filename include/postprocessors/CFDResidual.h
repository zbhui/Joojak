
#pragma once
#include "GeneralPostprocessor.h"

//Forward Declarations
class CFDResidual;

template<>
InputParameters validParams<CFDResidual>();

/**
 * Just returns the total number of Residual Evaluations performed.
 */
class CFDResidual : public GeneralPostprocessor
{
public:
	CFDResidual(const std::string & name, InputParameters parameters);

  virtual void initialize() {}
  virtual void execute() {}

  virtual Real getValue();
};
