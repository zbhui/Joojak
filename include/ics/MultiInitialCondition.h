
#pragma once

#include "InitialCondition.h"

class MultiInitialCondition :
	public InitialCondition
{
public:
	MultiInitialCondition(const std::string & name, InputParameters parameters);

	virtual Real value(const Point & p);
	virtual Real value(int component, const Point & p) = 0;

protected:
	unsigned int _component;
};

template<>
InputParameters validParams<MultiInitialCondition>();
