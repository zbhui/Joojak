
#include "CLawBoundaryCondition.h"

template<>
InputParameters validParams<CLawBoundaryCondition>()
{
	InputParameters params = validParams<IntegratedBC>();
	return params;
}

CLawBoundaryCondition::CLawBoundaryCondition(const std::string & name, InputParameters parameters):
		IntegratedBC(name, parameters),
		CLawInterface(parameters),
		_flux(getMaterialProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable(getMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable")),
		_flux_jacobi_grad_variable(getMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable")),

		_lift(getMaterialProperty<std::vector<RealVectorValue> >("lift")),
		_lift_jacobi_variable(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("lift_jacobi_variable_ee")),
		_eq(equationIndex(_var.name())),
		_epsilon(1),
		_sigma(6)
{
}

Real CLawBoundaryCondition::computeQpResidual()
{
	Real CIP = computeCIP();
	Real flux = _flux[_qp][_eq] ;
	flux += CIP*(_lift[_qp][_eq]+_lift[_qp][_eq])*_normals[_qp];
	return  flux * _test[_i][_qp] + _epsilon * _lift[_qp][_eq]* _grad_test[_i][_qp];
}

Real CLawBoundaryCondition::computeQpJacobian()
{
	Real r = 0;
	Real CIP = computeCIP();
	int p(_eq), q(_eq);
	r =  _flux_jacobi_variable[_qp][p][q]*_phi[_j][_qp]*_test[_i][_qp];
	r += _flux_jacobi_grad_variable[_qp][p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
	r += CIP*(_lift_jacobi_variable[_qp][p][q] + _lift_jacobi_variable[_qp][p][q])*_normals[_qp]*_phi[_j][_qp]*_test[_i][_qp];
	r += _epsilon*_lift_jacobi_variable[_qp][p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];

	return r;
}

Real CLawBoundaryCondition::computeQpOffDiagJacobian(unsigned int jvar)
{
	Real r = 0;
	Real CIP = computeCIP();
	int p(_eq), q(jvar);
	r =  _flux_jacobi_variable[_qp][p][q]*_phi[_j][_qp]*_test[_i][_qp];
	r += _flux_jacobi_grad_variable[_qp][p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
	r += CIP*(_lift_jacobi_variable[_qp][p][q] + _lift_jacobi_variable[_qp][p][q])*_normals[_qp]*_phi[_j][_qp]*_test[_i][_qp];
	r += _epsilon*_lift_jacobi_variable[_qp][p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];

	return r;
}

Real CLawBoundaryCondition::computeCIP()
{
	const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
	const double h_elem = (_current_elem_volume+_current_elem_volume)/_current_side_volume * 1./std::pow(elem_b_order, 2.)/2.;
	return _sigma/h_elem;
}
