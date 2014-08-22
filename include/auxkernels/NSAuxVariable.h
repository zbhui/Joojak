
#pragma once

#include "AuxKernel.h"
#include "CFDBase.h"

class NSAuxVariable;

template<>
InputParameters validParams<NSAuxVariable>();

/**
 * Coupled auxiliary value
 */
class NSAuxVariable :
public AuxKernel,
public CFDBase
{
public:

  /**
   * Factory constructor, takes parameters so that all derived classes can be built using the same
   * constructor.
   */
  NSAuxVariable(const std::string & name, InputParameters parameters);

protected:
    virtual Real computeValue();
	void valueAtCellPoint(Real *uh);

protected:
	std::vector<VariableValue*> _uh;
private:
};

