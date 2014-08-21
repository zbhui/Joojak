
#include "EulerBC.h"

template<>
InputParameters validParams<EulerBC>()
{
	InputParameters params = validParams<CFDBC>();

	return params;
}

EulerBC::EulerBC(const std::string & name, InputParameters parameters):
		CFDBC(name, parameters),
		_flux(getMaterialProperty<std::vector<Real> >("flux")),
		_jacobi_variable(getMaterialProperty<std::vector<std::vector<Real> > >("bnd_jacobi_variable"))
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

Real EulerBC::computeQpResidual()
{
	Real flux = _flux[_qp][_eq];
	return flux * _test[_i][_qp];
}

Real EulerBC::computeQpJacobian()
{
	return _jacobi_variable[_qp][_eq][_eq]*_phi[_j][_qp]*_test[_i][_qp];
}

Real EulerBC::computeQpOffDiagJacobian(unsigned int jvar)
{
	return _jacobi_variable[_qp][_eq][jvar]*_phi[_j][_qp]*_test[_i][_qp];
}
