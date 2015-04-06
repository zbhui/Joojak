
#pragma once

#include "MultiInitialCondition.h"

class PassFlowIC : public MultiInitialCondition
{
public:
	PassFlowIC(const std::string & name, InputParameters parameters);
	virtual Real value(int component, const Point & p);

};

template<>
InputParameters validParams<PassFlowIC>();
