
#pragma once

#include "ElementIntegralPostprocessor.h"
#include "EulerBase.h"

class BumpElementL2Error;

template<>
InputParameters validParams<BumpElementL2Error>();

class BumpElementL2Error :
public ElementIntegralPostprocessor,
public EulerBase
{
public:
	BumpElementL2Error(const std::string & name, InputParameters parameters);

	virtual Real getValue();

protected:
	virtual Real computeQpIntegral();
	PostprocessorValue &_area;
};
