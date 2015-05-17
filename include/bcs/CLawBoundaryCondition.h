
#pragma once

#include "IntegratedBC.h"
#include "CLawBoundaryMaterial.h"

class CLawBoundaryCondition : public IntegratedBC
{
public:
	CLawBoundaryCondition(const std::string & name, InputParameters params);

protected:
	virtual Real computeQpResidual();
	virtual Real computeQpJacobian();
	virtual Real computeQpOffDiagJacobian(unsigned int jvar);
	MaterialProperty<CLawBoundaryMaterialData> &_boundary;

private:
	Real computeQpJacobian(int p, int q);

	int _eq;
};

template<>
InputParameters validParams<CLawBoundaryCondition>();
