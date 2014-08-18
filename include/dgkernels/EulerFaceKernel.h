
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
	MaterialProperty<std::vector<RealVectorValue> > &_invis_term;
	MaterialProperty<std::vector<RealVectorValue> > &_invis_term_neighbor;
	MaterialProperty<Real > &_flux_diff;

	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);

	int _eq;
};
