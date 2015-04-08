
#pragma once

#include "FEProblem.h"

class CLawProblem : public FEProblem
{
public:
	CLawProblem(const std::string & name, InputParameters params);

	virtual int equationIndex(const std::string &var_name);

	virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh){};
	virtual void viscousTerm(RealVectorValue *viscous_term, Real* uh, RealGradient *duh){};
};

template<>
InputParameters validParams<CLawProblem>();
