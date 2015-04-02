

#pragma once

#include "InputParameters.h"

class CFDPostprocessorInterface
{
public:
	CFDPostprocessorInterface(InputParameters &parameter);
	~CFDPostprocessorInterface(){}

	void force();
	void momentun();

protected:
};
