
#pragma once

#include "IntegratedBC.h"
#include "EulerBase.h"

class EulerBC;

template<>
InputParameters validParams<EulerBC>();

class EulerBC :
public IntegratedBC,
public EulerBase
{
public:
	EulerBC(const std::string & name, InputParameters params);

protected:
	MaterialProperty<std::vector<Real> > &_flux;
	MaterialProperty<std::vector<std::vector<Real> > > &_jacobi_variable;
	virtual Real computeQpResidual();
	virtual Real computeQpJacobian();
	virtual Real computeQpOffDiagJacobian(unsigned int jvar);

	int _eq;
};
