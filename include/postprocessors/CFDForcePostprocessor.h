
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
	enum ForceType
	{
		form = 0,
		friction = 1,
		total = 2
	};
	enum Direction
	{
		x = 0,
		y = 1,
		z = 2
	};
	Direction _direction;
	ForceType _force_type;

	int _n_equations;
	std::vector<VariableValue*> _uh;
	std::vector<VariableGradient*> _grad_uh;
	virtual Real computeQpIntegral();
	void computeQpValue(Real *uh, RealGradient *duh);
};
