#pragma once

#include "Kernel.h"
#include "CLawInterface.h"
#include "CLawCellMaterial.h"

class CLawProblem;

class CLawCellKernel :
public Kernel,
public CLawInterface
{
public:
	CLawCellKernel(const std::string & name, InputParameters parameters);
	virtual ~CLawCellKernel(){}

protected:
	int _eq;
	MaterialProperty<CLawCellMaterialData > &_cell_material_data;

	virtual Real computeQpResidual();
	virtual Real computeQpJacobian();
	virtual Real computeQpOffDiagJacobian(unsigned int jvar);

private:
	Real computeQpJacobian(int p, int q);
};

template<>
InputParameters validParams<CLawCellKernel>();
