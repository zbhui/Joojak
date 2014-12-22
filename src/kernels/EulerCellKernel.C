#include "EulerCellKernel.h"

template<>
InputParameters validParams<EulerCellKernel>()
{
  InputParameters params = validParams<Kernel>();
  params += validParams<EulerBase>();
  return params;
}

EulerCellKernel::EulerCellKernel(const std::string & name, InputParameters parameters):
		Kernel(name, parameters),
		EulerBase(name, parameters),
		_invis_term(getMaterialProperty<std::vector<RealVectorValue> >("cell_material")),
		_jacobi(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("cell_jacobi_variable"))
{
	_eq = equationIndex(_var.name());
}

Real EulerCellKernel::computeQpResidual()
{
	return -(_invis_term[_qp][_eq]*_grad_test[_i][_qp]);
}

Real EulerCellKernel::computeQpJacobian()
{
	return -_jacobi[_qp][_eq][_eq]*_phi[_j][_qp]*_grad_test[_i][_qp];
}
Real EulerCellKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
	return -_jacobi[_qp][_eq][jvar]*_phi[_j][_qp]*_grad_test[_i][_qp];
}
