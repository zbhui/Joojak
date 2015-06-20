
#pragma once

#include "DGKernel.h"
#include "CLawFaceMaterial.h"

class CLawFaceKernel :
public DGKernel
{
public:
	CLawFaceKernel(const std::string &name, InputParameters parameters);

protected:
	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);
	virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar);

	int _eq;
	const MaterialProperty<CLawFaceMaterialData> &_face;

private:
	Real computeQpJacobian(int p, int q, Moose::DGJacobianType type);
};

template<>
InputParameters validParams<CLawFaceKernel>();
