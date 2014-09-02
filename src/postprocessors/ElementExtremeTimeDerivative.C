
#include "ElementExtremeTimeDerivative.h"

template<>
InputParameters validParams<ElementExtremeTimeDerivative>()
{
  InputParameters params = validParams<ElementExtremeValue>();
  return params;
}

ElementExtremeTimeDerivative::ElementExtremeTimeDerivative(const std::string & name, InputParameters parameters) :
		ElementExtremeValue(name, parameters)
{}

void
ElementExtremeTimeDerivative::computeQpValue()
{
  switch (_type)
  {
    case MAX:
      _value = std::max(_value, std::fabs(_u_dot[_qp]));
      break;

    case MIN:
      _value = std::min(_value, std::fabs(_u_dot[_qp]));
      break;
  }
}

