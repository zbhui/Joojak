
#pragma once

#include "CFDBC.h"
#include "NSBase.h"

class NSBC;

template<>
InputParameters validParams<NSBC>();

class NSBC :
public CFDBC,
public NSBase
{
public:
	NSBC(const std::string & name, InputParameters params);

protected:
	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > &_flux_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealGradient> > > &_flux_jacobi_grad_variable;
	MaterialProperty<std::vector<RealVectorValue> > & _penalty;
	MaterialProperty<std::vector<RealVectorValue> > & _penalty_neighbor;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > &_penalty_jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > &_penalty_jacobi_variable_ne;

	virtual Real computeQpResidual();
	virtual Real computeQpJacobian();
	virtual Real computeQpOffDiagJacobian(unsigned int jvar);

	Real computeCIP();
	int _eq;
};
