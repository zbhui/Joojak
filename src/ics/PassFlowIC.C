

#include "PassFlowIC.h"

template<>
InputParameters validParams<PassFlowIC>()
{
	InputParameters params = validParams<MultiInitialCondition>();
	return params;
}

PassFlowIC::PassFlowIC(const std::string& name, InputParameters parameters) :
		MultiInitialCondition(name, parameters)
{
//	_fe_problem
}

Real PassFlowIC::value(int component, const Point& p)
{
	return 0;
}
