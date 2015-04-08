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
//		_claw_problem(static_cast<CLawProblem&>(_fe_problem)),
		_eq(_claw_problem.equationIndex((_var.name()))),
		_flux_term(getMaterialProperty<std::vector<RealVectorValue> >("flux_term")),
		_flux_jacobi_variable(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("flux_term_jacobi_variable")),
		_flux_jacobi_grad_variable(getMaterialProperty<std::vector<std::vector<RealTensorValue> > >("flux_term_jacobi_grad_variable"))
{
//	_eq = equationIndex(_var.name());
}

Real CLawCellKernel::computeQpResidual()
{
	return -(_flux_term[_qp][_eq]*_grad_test[_i][_qp]);
}

Real CLawCellKernel::computeQpJacobian()
{
	return -(_flux_jacobi_variable[_qp][_eq][_eq]*_phi[_j][_qp]+_flux_jacobi_grad_variable[_qp][_eq][_eq]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
}
Real CLawCellKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
	return -(_flux_jacobi_variable[_qp][_eq][jvar]*_phi[_j][_qp]+_flux_jacobi_grad_variable[_qp][_eq][jvar]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
}
