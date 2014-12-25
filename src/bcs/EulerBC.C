
#include "EulerBC.h"

template<>
InputParameters validParams<EulerBC>()
{
	InputParameters params = validParams<IntegratedBC>();
	params += validParams<EulerBase>();

	return params;
}

EulerBC::EulerBC(const std::string & name, InputParameters parameters):
		IntegratedBC(name, parameters),
		EulerBase(name, parameters),
		_flux(getMaterialProperty<std::vector<Real> >("flux")),
		_jacobi_variable(getMaterialProperty<std::vector<std::vector<Real> > >("bnd_jacobi_variable"))
{
	_eq = equationIndex(_var.name());
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
