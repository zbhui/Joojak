
#pragma once

#include "DGKernel.h"
#include "EulerBase.h"

class EulerFaceKernel;

template<>
InputParameters validParams<EulerFaceKernel>();

class EulerFaceKernel :
public DGKernel,
public EulerBase
{
public:
	EulerFaceKernel(const std::string &name, InputParameters parameters);
	virtual ~EulerFaceKernel(){}

protected:
	MaterialProperty<std::vector<Real> > &_flux;
	MaterialProperty<std::vector<std::vector<Real> > > & _jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<Real> > > & _jacobi_variable_en;
	MaterialProperty<std::vector<std::vector<Real> > > & _jacobi_variable_ne;
	MaterialProperty<std::vector<std::vector<Real> > > & _jacobi_variable_nn;

	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);
	virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar);

	int _eq;
};
