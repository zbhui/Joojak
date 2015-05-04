
#pragma once

#include "ElementIntegralPostprocessor.h"
#include "IsoVortexBase.h"

class IsoVortexElementL2Error :
public ElementIntegralPostprocessor,
public IsoVortexBase
{
public:
	IsoVortexElementL2Error(const std::string & name, InputParameters parameters);

	virtual Real getValue();

protected:
	virtual Real computeQpIntegral();

private:
	std::vector<VariableValue*> _uh;
};

template<>
InputParameters validParams<IsoVortexElementL2Error>();
