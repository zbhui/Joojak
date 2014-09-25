#include "SACellKernel.h"

template<>
InputParameters validParams<SACellKernel>()
{
  InputParameters params = validParams<Kernel>();
  params += validParams<SABase>();
  return params;
}

SACellKernel::SACellKernel(const std::string & name, InputParameters parameters):
		Kernel(name, parameters),
		SABase(name, parameters),
		_flux_term(getMaterialProperty<std::vector<RealVectorValue> >("flux_term")),
		_source_term(getMaterialProperty<std::vector<Real> >("source_term")),
		_flux_jacobi_variable(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("flux_term_jacobi_variable")),
		_flux_jacobi_grad_variable(getMaterialProperty<std::vector<std::vector<RealTensorValue> > >("flux_term_jacobi_grad_variable")),
		_source_jacobi_variable(getMaterialProperty<std::vector<std::vector<Real> > >("source_term_jacobi_variable")),
		_source_jacobi_grad_variable(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("source_term_jacobi_grad_variable"))
{
	_eq = equationIndex(_var.name());
}

Real SACellKernel::computeQpResidual()
{
	return -(_flux_term[_qp][_eq]*_grad_test[_i][_qp])-_source_term[_qp][_eq]*_test[_i][_qp];
}

Real SACellKernel::computeQpJacobian()
{
	int p(_eq), q(_eq);
	Real r = -(_flux_jacobi_variable[_qp][p][q]*_phi[_j][_qp]+_flux_jacobi_grad_variable[_qp][p][q]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
	r += -(_source_jacobi_variable[_qp][p][q]*_phi[_j][_qp]+_source_jacobi_grad_variable[_qp][p][q]*_grad_phi[_j][_qp])*_test[_i][_qp];
	return r;
}

Real SACellKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
	int p(_eq), q(jvar);
	Real r = -(_flux_jacobi_variable[_qp][p][q]*_phi[_j][_qp]+_flux_jacobi_grad_variable[_qp][p][q]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
	r += -(_source_jacobi_variable[_qp][p][q]*_phi[_j][_qp]+_source_jacobi_grad_variable[_qp][p][q]*_grad_phi[_j][_qp])*_test[_i][_qp];
	return r;
}
