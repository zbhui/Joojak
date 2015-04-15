#include "CLawCellKernel.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawCellKernel>()
{
  InputParameters params = validParams<Kernel>();
  return params;
}

CLawCellKernel::CLawCellKernel(const std::string & name, InputParameters parameters):
		Kernel(name, parameters),
		CLawInterface(parameters),
		_eq(equationIndex((_var.name()))),
		_flux_term(getMaterialProperty<std::vector<RealVectorValue> >("flux_term")),
		_flux_jacobi_variable(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("flux_term_jacobi_variable")),
		_flux_jacobi_grad_variable(getMaterialProperty<std::vector<std::vector<RealTensorValue> > >("flux_term_jacobi_grad_variable")),
		_cell_material_data(getMaterialProperty<CLawCellMaterialData>("cell_material_data"))
{
}

Real CLawCellKernel::computeQpResidual()
{
	RealVectorValue flux = _cell_material_data[_qp]._flux_term[_eq];
	return -(flux*_grad_test[_i][_qp]);
//	return -(_flux_term[_qp][_eq]*_grad_test[_i][_qp]);
}

Real CLawCellKernel::computeQpJacobian()
{
	return -computeQpJacobian(_eq, _eq);
}
Real CLawCellKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
	return -computeQpJacobian(_eq, jvar);
}

Real CLawCellKernel::computeQpJacobian(int p, int q)
{
	return (_flux_jacobi_variable[_qp][p][q]*_phi[_j][_qp]+_flux_jacobi_grad_variable[_qp][p][q]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];

}
