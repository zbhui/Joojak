

#pragma once

#include "InputParameters.h"

class CLawProblem;

class CLawInterface
{
public:
	CLawInterface(InputParameters &parameter);
	~CLawInterface(){}

	void convertionTerm();
	void diffusionTerm();
	void sourceTerm();

protected:
	CLawProblem &_claw_problem;
};
