
#include "NSBC.h"

template<>
InputParameters validParams<NSBC>()
{
	InputParameters params = validParams<CFDBC>();
	params += validParams<NSBase>();

	return params;
}

NSBC::NSBC(const std::string & name, InputParameters parameters):
		CFDBC(name, parameters),
		NSBase(name, parameters),
		_flux(getMaterialProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable(getMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable")),
		_flux_jacobi_grad_variable(getMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable")),

		_penalty(getMaterialProperty<std::vector<RealVectorValue> >("penalty")),
		_penalty_jacobi_variable(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable"))
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

Real NSBC::computeQpResidual()
{
	Real CIP = computeCIP();
	Real flux = _flux[_qp][_eq] ;
	flux += CIP*(_penalty[_qp][_eq] + _penalty[_qp][_eq])*_normals[_qp];
	return  _flux[_qp][_eq] * _test[_i][_qp] + 0.5*_epsilon * _penalty[_qp][_eq]* _grad_test[_i][_qp];
}

Real NSBC::computeQpJacobian()
{
	return 0.;
//	return _flux_jacobi_variable[_qp][_eq][_eq]*_phi[_j][_qp]*_test[_i][_qp];
}

Real NSBC::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.;
//	return _flux_jacobi_variable[_qp][_eq][jvar]*_phi[_j][_qp]*_test[_i][_qp];
}

Real NSBC::computeCIP()
{
	const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
	const double h_elem = (_current_elem_volume+_current_elem_volume)/_current_side_volume * 1./std::pow(elem_b_order, 2.)/2.;
	return _sigma/h_elem;
}
