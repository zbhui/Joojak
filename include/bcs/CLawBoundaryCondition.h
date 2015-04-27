
#pragma once

#include "IntegratedBC.h"
#include "CLawInterface.h"

class CLawBoundaryCondition :
public IntegratedBC,
public CLawInterface
{
public:
	CLawBoundaryCondition(const std::string & name, InputParameters params);

protected:
	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > &_flux_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealGradient> > > &_flux_jacobi_grad_variable;
	MaterialProperty<std::vector<RealVectorValue> > & _lift;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > &_lift_jacobi_variable;

	virtual Real computeQpResidual();
	virtual Real computeQpJacobian();
	virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
	Real computeQpJacobian(int p, int q);

	int _eq;
};

template<>
InputParameters validParams<CLawBoundaryCondition>();
