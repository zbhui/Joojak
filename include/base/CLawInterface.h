

#pragma once

#include "InputParameters.h"

class CLawProblem;

class CLawInterface
{
public:
	CLawInterface(InputParameters &parameter);
	~CLawInterface(){}

	void convertionTerm(RealVectorValue *inviscous_term, Real *uh);
	void diffusionTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh);
	void sourceTerm();

protected:
	CLawProblem &_claw_problem;
};
