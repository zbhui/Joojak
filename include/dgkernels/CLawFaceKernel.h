
#pragma once

#include "DGKernel.h"
#include "CLawInterface.h"

class CLawFaceKernel :
public DGKernel,
public CLawInterface
{
public:
	CLawFaceKernel(const std::string &name, InputParameters parameters);

protected:
	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_en;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_ee;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_en;

	MaterialProperty<std::vector<RealVectorValue> > & _penalty;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _penalty_jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _penalty_jacobi_variable_en;

	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);
	virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar);

	int _eq;
	Real _epsilon;
};

template<>
InputParameters validParams<CLawFaceKernel>();
