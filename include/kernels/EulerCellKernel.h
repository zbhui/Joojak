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
	int _eq;
	MaterialProperty<std::vector<RealVectorValue> > &_invis_term;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > >& _jacobi;

	virtual Real computeQpResidual();
	virtual Real computeQpJacobian();
	virtual Real computeQpOffDiagJacobian(unsigned int jvar);

};
