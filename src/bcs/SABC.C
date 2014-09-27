
#include "SABC.h"

template<>
InputParameters validParams<SABC>()
{
	MooseEnum bc_types("wall far_field symmetric pressure_out none", "none");  // 边界条件的类型，可以增加

	InputParameters params = validParams<IntegratedBC>();
	params += validParams<SABase>();
	params.addRequiredParam<MooseEnum>("bc_type", bc_types, "边界条件");
	return params;
}

SABC::SABC(const std::string & name, InputParameters parameters):
		IntegratedBC(name, parameters),
		SABase(name, parameters),
		_bc_type(getParam<MooseEnum>("bc_type")),
		_flux(getMaterialProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable(getMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable")),
		_flux_jacobi_grad_variable(getMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable")),

		_penalty(getMaterialProperty<std::vector<RealVectorValue> >("penalty")),
		_penalty_neighbor(getMaterialProperty<std::vector<RealVectorValue> >("penalty_neighbor")),
		_penalty_jacobi_variable_ee(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_ee")),
		_penalty_jacobi_variable_ne(getMaterialProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_ne"))
{
	_eq = equationIndex(_var.name());
}

Real SABC::computeQpResidual()
{
	Real CIP = computeCIP();
	Real flux = _flux[_qp][_eq] ;
	flux += CIP*(_penalty[_qp][_eq] + _penalty_neighbor[_qp][_eq])*_normals[_qp];
	return  flux * _test[_i][_qp] + _epsilon * _penalty[_qp][_eq]* _grad_test[_i][_qp];
}

Real SABC::computeQpJacobian()
{
	Real r = 0;
	Real CIP = computeCIP();
	int p(_eq), q(_eq);
	r =  _flux_jacobi_variable[_qp][p][q]*_phi[_j][_qp]*_test[_i][_qp];
	r += _flux_jacobi_grad_variable[_qp][p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
	r += CIP*(_penalty_jacobi_variable_ee[_qp][p][q] + _penalty_jacobi_variable_ne[_qp][p][q])*_normals[_qp]*_phi[_j][_qp]*_test[_i][_qp];
	r += _epsilon*_penalty_jacobi_variable_ee[_qp][p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];

	return r;
}

Real SABC::computeQpOffDiagJacobian(unsigned int jvar)
{
	Real r = 0;
	Real CIP = computeCIP();
	int p(_eq), q(jvar);
	r =  _flux_jacobi_variable[_qp][p][q]*_phi[_j][_qp]*_test[_i][_qp];
	r += _flux_jacobi_grad_variable[_qp][p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
	r += CIP*(_penalty_jacobi_variable_ee[_qp][p][q] + _penalty_jacobi_variable_ne[_qp][p][q])*_normals[_qp]*_phi[_j][_qp]*_test[_i][_qp];
	r += _epsilon*_penalty_jacobi_variable_ee[_qp][p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];
	return r;
}

Real SABC::computeCIP()
{
	const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
	const double h_elem = (_current_elem_volume+_current_elem_volume)/_current_side_volume * 1./std::pow(elem_b_order, 2.)/2.;
	return _sigma/h_elem;
}
