#include "EulerCellKernel.h"

template<>
InputParameters validParams<EulerCellKernel>()
{
  InputParameters params = validParams<Kernel>();
//  params.addRequiredParam<std::string>("material", "单元的材料属性");
  return params;
}

EulerCellKernel::EulerCellKernel(const std::string & name, InputParameters parameters):
		Kernel(name, parameters),
		_invis_term(getMaterialProperty<std::vector<RealVectorValue> >("cell_material")),
		_jacobi(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("cell_jacobi_variable"))
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
