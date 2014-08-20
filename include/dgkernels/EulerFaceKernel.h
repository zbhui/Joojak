
#pragma once

#include "DGKernel.h"

class EulerFaceKernel;

template<>
InputParameters validParams<EulerFaceKernel>();

class EulerFaceKernel :
public DGKernel
{
public:
	EulerFaceKernel(const std::string &name, InputParameters parameters);
	virtual ~EulerFaceKernel(){}

protected:
	MaterialProperty<std::vector<Real> > &_flux;

	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);

	int _eq;
};
