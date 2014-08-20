#pragma once

#include "Kernel.h"

class EulerCellKernel;

template<>
InputParameters validParams<EulerCellKernel>();

class EulerCellKernel :
public Kernel
{
public:
	EulerCellKernel(const std::string & name, InputParameters parameters);
	virtual ~EulerCellKernel(){}

protected:
	MaterialProperty<std::vector<RealVectorValue> > &_invis_term;

	virtual Real computeQpResidual();
	virtual Real computeQpJacobian();

	int _eq;
};
