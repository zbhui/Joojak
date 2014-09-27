#pragma once

#include "Kernel.h"
#include "SABase.h"

class SACellKernel;

template<>
InputParameters validParams<SACellKernel>();

class SACellKernel :
public Kernel,
public SABase
{
public:
	SACellKernel(const std::string & name, InputParameters parameters);

protected:
	int _eq;
	MaterialProperty<std::vector<RealVectorValue> > &_flux_term;
	MaterialProperty<std::vector<Real> > & _source_term;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > >& _flux_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealTensorValue> > >& _flux_jacobi_grad_variable;
	MaterialProperty<std::vector<std::vector<Real> > >& _source_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > >& _source_jacobi_grad_variable;


	virtual Real computeQpResidual();
	virtual Real computeQpJacobian();
	virtual Real computeQpOffDiagJacobian(unsigned int jvar);

};
