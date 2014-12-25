#include "NSCellKernel.h"

template<>
InputParameters validParams<NSCellKernel>()
{
  InputParameters params = validParams<Kernel>();
  params += validParams<NSBase>();
  return params;
}

NSCellKernel::NSCellKernel(const std::string & name, InputParameters parameters):
		Kernel(name, parameters),
		NSBase(name, parameters),
		_flux_term(getMaterialProperty<std::vector<RealVectorValue> >("flux_term")),
		_flux_jacobi_variable(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("flux_term_jacobi_variable")),
		_flux_jacobi_grad_variable(getMaterialProperty<std::vector<std::vector<RealTensorValue> > >("flux_term_jacobi_grad_variable"))
{
	_eq = equationIndex(_var.name());
}

Real NSCellKernel::computeQpResidual()
{
	return -(_flux_term[_qp][_eq]*_grad_test[_i][_qp]);
}

Real NSCellKernel::computeQpJacobian()
{
	return -(_flux_jacobi_variable[_qp][_eq][_eq]*_phi[_j][_qp]+_flux_jacobi_grad_variable[_qp][_eq][_eq]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
}
Real NSCellKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
	return -(_flux_jacobi_variable[_qp][_eq][jvar]*_phi[_j][_qp]+_flux_jacobi_grad_variable[_qp][_eq][jvar]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
}
