
#pragma once

#include "InputParameters.h"
#include "MooseTypes.h"

#include <string>

class MooseMesh;
class SubProblem;

class CLawBase
{
public:
	CLawBase(const std::string &name, InputParameters &paramter);
	virtual ~CLawBase(){}

	virtual void convertionTerm() = 0;
	virtual void diffusionTerm() = 0;
	virtual void sourceTerm() = 0;

protected:
	std::string _name;
	SubProblem & _subproblem;
	MooseMesh & _mesh;
};
