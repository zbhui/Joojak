
#pragma once

#include "AuxKernel.h"
#include "NSBase.h"

class NSAuxVariable;

template<>
InputParameters validParams<NSAuxVariable>();

class NSAuxVariable :
public AuxKernel,
public NSBase
{
public:
  NSAuxVariable(const std::string & name, InputParameters parameters);

protected:
  virtual Real computeValue();
  void valueAtCellPoint(Real *uh);
};

