
#pragma once

#include "FEProblem.h"

class CLawProblem : public FEProblem
{
public:
	CLawProblem(const std::string & name, InputParameters params);

	virtual int equationIndex(const std::string &var_name);
};

template<>
InputParameters validParams<CLawProblem>();
