

#include "MultiInitialCondition.h"

template<>
InputParameters validParams<MultiInitialCondition>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<int>("component", "分量");
  return params;
}

MultiInitialCondition::MultiInitialCondition(const std::string & name, InputParameters parameters) :
    InitialCondition(name, parameters),
	_component(getParam<int>("component"))
{}

Real MultiInitialCondition::value(const Point & p)
{
  return value(_component, p);
}
