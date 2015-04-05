
#pragma once

#include "FEProblem.h"

class CLawProblem : public FEProblem
{
public:
	CLawProblem(const std::string & name, InputParameters params);
};

template<>
InputParameters validParams<CLawProblem>();
