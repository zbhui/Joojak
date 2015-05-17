#include "CLawCellKernel.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawCellKernel>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<int>("component", "Kernel component");
  return params;
}

CLawCellKernel::CLawCellKernel(const std::string & name, InputParameters parameters):
		Kernel(name, parameters),
		_eq(getParam<int>("component")),
		_cell(getMaterialProperty<CLawCellMaterialData>("cell_material_data"))
{
}

Real CLawCellKernel::computeQpResidual()
{
	return -(_cell[_qp]._flux_term[_eq]*_grad_test[_i][_qp]) -_cell[_qp]._source_term[_eq]*_test[_i][_qp];
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
	Real r = (_cell[_qp]._flux_jacobi_variable[p][q]*_phi[_j][_qp]+_cell[_qp]._flux_jacobi_grad_variable[p][q]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
	r += (_cell[_qp]._source_jacobi_variable[p][q]*_phi[_j][_qp] + _cell[_qp]._source_jacobi_grad_variable[p][q]*_grad_phi[_j][_qp])*_test[_i][_qp];
	return r ;
}
