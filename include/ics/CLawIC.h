
#pragma once

#include "MultiInitialCondition.h"
class CLawProblem;

class CLawIC : public MultiInitialCondition
{
public:
	CLawIC(const std::string & name, InputParameters parameters);

	virtual Real value(int component, const Point & p);

private:
	CLawProblem &_claw_problem;
};

template<>
InputParameters validParams<CLawIC>();
