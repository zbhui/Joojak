
#pragma once

#include "DGKernel.h"
#include "EulerMaterial.h"

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
	MaterialProperty<RealVectorValue*> &_invis_term;
	MaterialProperty<RealVectorValue*> &_invis_term_neighbor;

	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);

	int _eq;
};
