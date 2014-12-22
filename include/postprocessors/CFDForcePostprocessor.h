
#pragma once

#include "SideIntegralPostprocessor.h"
#include "NSBase.h"

class CFDForcePostprocessor;

template<>
InputParameters validParams<CFDForcePostprocessor>();

class CFDForcePostprocessor :
public SideIntegralPostprocessor,
public NSBase
{
public:
	CFDForcePostprocessor(const std::string & name, InputParameters parameters);
	virtual void threadJoin(const UserObject &y);

protected:
	MooseEnum _direction;
	MooseEnum _force_type;

	int _n_equations;
	std::vector<VariableValue*> _uh;
	std::vector<VariableGradient*> _grad_uh;
	virtual Real computeQpIntegral();
	void computeQpValue(Real *uh, RealGradient *duh);
};
