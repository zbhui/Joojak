
#pragma once

#include "ElementExtremeValue.h"

class ElementExtremeTimeDerivative;

template<>
InputParameters validParams<ElementExtremeTimeDerivative>();

class ElementExtremeTimeDerivative :
public ElementExtremeValue
{
public:
	ElementExtremeTimeDerivative(const std::string & name, InputParameters parameters);

protected:
  virtual void computeQpValue();
};
