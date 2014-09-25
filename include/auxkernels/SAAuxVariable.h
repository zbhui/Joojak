
#pragma once

#include "AuxKernel.h"
#include "SABase.h"

class SAAuxVariable;

template<>
InputParameters validParams<SAAuxVariable>();

/**
 * Coupled auxiliary value
 */
class SAAuxVariable :
public AuxKernel,
public SABase
{
public:
	SAAuxVariable(const std::string & name, InputParameters parameters);

protected:
    virtual Real computeValue();
	void valueAtCellPoint(Real *uh);
};

