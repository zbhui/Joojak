

#pragma once

#include "InputParameters.h"

class CLawBase;

class CLawInterface
{
public:
	CLawInterface(InputParameters &parameter);
	~CLawInterface(){}

	void convertionTerm();
	void diffusionTerm();
	void sourceTerm();

protected:
	CLawBase &_claw_base;
};
