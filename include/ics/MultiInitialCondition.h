
#pragma once

#include "InitialCondition.h"

typedef Real (*buildICs)(const Point & p);

class MultiInitialCondition :
	public InitialCondition
{
public:
	MultiInitialCondition(const std::string & name, InputParameters parameters);

	virtual Real value(const Point & p);
	virtual Real value(int component, const Point & p) = 0;

protected:
	unsigned int _component;
	std::map<std::string, buildICs> _name_to_ics;
};

template<>
InputParameters validParams<MultiInitialCondition>();
