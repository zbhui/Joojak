
#pragma once

#include "IntegratedBC.h"
#include "MooseEnum.h"
#include "SABase.h"

class SABC;

template<>
InputParameters validParams<SABC>();

class SABC :
public IntegratedBC,
public SABase
{
public:
	SABC(const std::string & name, InputParameters params);

protected:
	MooseEnum _bc_type;
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
