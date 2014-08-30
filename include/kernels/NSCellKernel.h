#pragma once

#include "Kernel.h"

class NSCellKernel;

template<>
InputParameters validParams<NSCellKernel>();

class NSCellKernel :
public Kernel
{
public:
	NSCellKernel(const std::string & name, InputParameters parameters);

protected:
	int _eq;
	MaterialProperty<std::vector<RealVectorValue> > &_flux_term;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > >& _flux_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealTensorValue> > >& _flux_jacobi_grad_variable;

	virtual Real computeQpResidual();
	virtual Real computeQpJacobian();
	virtual Real computeQpOffDiagJacobian(unsigned int jvar);

};
