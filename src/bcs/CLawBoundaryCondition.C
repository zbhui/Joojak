
#include "CLawBoundaryCondition.h"

template<>
InputParameters validParams<CLawBoundaryCondition>()
{
	InputParameters params = validParams<IntegratedBC>();
	params.addRequiredParam<int>("component", "BC component");
	return params;
}

CLawBoundaryCondition::CLawBoundaryCondition(const std::string & name, InputParameters parameters):
		IntegratedBC(name, parameters),
//		CLawInterface(parameters),
		_flux(getMaterialProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable(getMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable")),
		_flux_jacobi_grad_variable(getMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable")),

		_lift(getMaterialProperty<std::vector<RealVectorValue> >("lift")),
		_lift_jacobi_variable(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("lift_jacobi_variable")),
		_eq(getParam<int>("component"))
{
}

Real CLawBoundaryCondition::computeQpResidual()
{
	return _flux[_qp][_eq] * _test[_i][_qp] + _lift[_qp][_eq]*_grad_test[_i][_qp];
}

Real CLawBoundaryCondition::computeQpJacobian()
{
	int p(_eq), q(_eq);
	return computeQpJacobian(p, q);
}

Real CLawBoundaryCondition::computeQpOffDiagJacobian(unsigned int jvar)
{
	int p(_eq), q(jvar);
	return computeQpJacobian(p, q);
}

Real CLawBoundaryCondition::computeQpJacobian(int p, int q)
{
	Real r(0);

	r = _flux_jacobi_variable[_qp][p][q]*_phi[_j][_qp]*_test[_i][_qp];
	r += _flux_jacobi_grad_variable[_qp][p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
	r += _lift_jacobi_variable[_qp][p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];

	return r;
}
