
#include "EulerBC.h"

template<>
InputParameters validParams<EulerBC>()
{
	InputParameters params = validParams<CFDBC>();

	return params;
}

EulerBC::EulerBC(const std::string & name, InputParameters parameters):
		CFDBC(name, parameters),
		_invis_term(getMaterialProperty<std::vector<RealVectorValue> >("left_material")),
		_invis_term_neighbor(getMaterialProperty<std::vector<RealVectorValue> >("right_material")),
		_flux_diff(getMaterialProperty<Real>("flux_diff")),
		_uh_neighbor(getMaterialProperty<std::vector<Real> >("right_value"))
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
	Real flux = 0.5*(_invis_term[_qp][_eq] + _invis_term_neighbor[_qp][_eq])*_normals[_qp];
	flux += 1*(_u[_qp]-_uh_neighbor[_qp][_eq]);

	return flux * _test[_i][_qp];
}
