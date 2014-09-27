
#pragma once

#include "DGKernel.h"
#include "KOBase.h"

class KOFaceKernel;

template<>
InputParameters validParams<KOFaceKernel>();

class KOFaceKernel :
public DGKernel,
public KOBase
{
public:
	KOFaceKernel(const std::string &name, InputParameters parameters);

protected:
	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_en;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_ne;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_nn;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_ee;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_en;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_ne;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_nn;

	MaterialProperty<std::vector<RealVectorValue> > & _penalty;
	MaterialProperty<std::vector<RealVectorValue> > & _penalty_neighbor;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _penalty_jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _penalty_jacobi_variable_en;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _penalty_jacobi_variable_ne;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _penalty_jacobi_variable_nn;

	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);
	virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar);

	Real computeCIP();

	int _eq;
};
