
#pragma once

#include "DGKernel.h"

class CLawFaceKernel :
public DGKernel
{
public:
	CLawFaceKernel(const std::string &name, InputParameters parameters);

protected:
	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_en;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_ee;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_en;

	MaterialProperty<std::vector<RealVectorValue> > & _lift;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _lift_jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _lift_jacobi_variable_en;

	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);
	virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar);

	int _eq;

private:
	Real computeQpJacobian(int p, int q, Moose::DGJacobianType type);
};

template<>
InputParameters validParams<CLawFaceKernel>();
