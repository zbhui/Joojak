#include "KOCellKernel.h"

template<>
InputParameters validParams<KOCellKernel>()
{
  InputParameters params = validParams<Kernel>();
  return params;
}

KOCellKernel::KOCellKernel(const std::string & name, InputParameters parameters):
		Kernel(name, parameters),
		_flux_term(getMaterialProperty<std::vector<RealVectorValue> >("flux_term")),
		_flux_jacobi_variable(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("flux_term_jacobi_variable")),
		_flux_jacobi_grad_variable(getMaterialProperty<std::vector<std::vector<RealTensorValue> > >("flux_term_jacobi_grad_variable"))
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
	if(var_name == "rhok")
		_eq = 5;
	if(var_name == "rhoo")
		_eq = 6;
}

Real KOCellKernel::computeQpResidual()
{
	return -(_flux_term[_qp][_eq]*_grad_test[_i][_qp]);
}

Real KOCellKernel::computeQpJacobian()
{
	return -(_flux_jacobi_variable[_qp][_eq][_eq]*_phi[_j][_qp]+_flux_jacobi_grad_variable[_qp][_eq][_eq]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
}
Real KOCellKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
	return -(_flux_jacobi_variable[_qp][_eq][jvar]*_phi[_j][_qp]+_flux_jacobi_grad_variable[_qp][_eq][jvar]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
}
