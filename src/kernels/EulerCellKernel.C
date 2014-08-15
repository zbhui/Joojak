#include "EulerCellKernel.h"

template<>
InputParameters validParams<EulerCellKernel>()
{
  InputParameters params = validParams<Kernel>();

  return params;
}

EulerCellKernel::EulerCellKernel(const std::string & name, InputParameters parameters):
		Kernel(name, parameters),
		_inviscous_term(getMaterialProperty<RealVectorValue*>("inviscous"))
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

Real EulerCellKernel::computeQpJacobian()
{
	return 0.;
}

Real EulerCellKernel::computeQpResidual()
{
	return -_inviscous_term[_qp][_eq]*_grad_test[_i][_qp];
}
