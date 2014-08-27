#include "NSCellKernel.h"

template<>
InputParameters validParams<NSCellKernel>()
{
  InputParameters params = validParams<Kernel>();
  return params;
}

NSCellKernel::NSCellKernel(const std::string & name, InputParameters parameters):
		Kernel(name, parameters),
		_flux_term(getMaterialProperty<std::vector<RealVectorValue> >("flux_term")),
		_jacobi_variable(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("jacobi_variable")),
		_jacobi_grad_variable(getMaterialProperty<std::vector<std::vector<RealTensorValue> > >("jacobi_grad_variable"))
{
	std::string var_name = _var.name();

	if(var_name == "rho")
		_eq = 0;
	if(var_name == "momentum_x")
		_eq = 1;
	if(var_name == "momentum_y")
		_eq = 2;
	if(var_name == "momentum_z")
		_eq = 3;
	if(var_name == "rhoe")
		_eq = 4;
}

Real NSCellKernel::computeQpResidual()
{
	return -(_flux_term[_qp][_eq]*_grad_test[_i][_qp]);
}

Real NSCellKernel::computeQpJacobian()
{
	return 0.;
//	return -_jacobi[_qp][_eq][_eq]*_phi[_j][_qp]*_grad_test[_i][_qp];
}
Real NSCellKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.;
//	return -_jacobi[_qp][_eq][jvar]*_phi[_j][_qp]*_grad_test[_i][_qp];
}
