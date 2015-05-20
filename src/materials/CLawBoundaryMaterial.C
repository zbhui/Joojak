
#include "CLawBoundaryMaterial.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawBoundaryMaterial>()
{
  InputParameters params = validParams<CLawMaterial>();

  params.addParam<std::string>("bc_type", "empty", "边界条件");
  params.addParam<Real>("ds", 1.490116119384766e-08, "微扰量");
  params.addParam<Real>("sigma", 6, "通量罚值，默认值为6");
  params.addParam<Real>("epsilon", 1, "对称项罚值，可以取1, 0 , -1，分别对应SIP, IIP, NIP");

  return params;
}

CLawBoundaryMaterial::CLawBoundaryMaterial(const std::string & name, InputParameters parameter):
		CLawMaterial(name, parameter),
		_bc_type(getParam<std::string>("bc_type")),
		_ds(getParam<Real>("ds")),
		_sigma(getParam<Real>("sigma")),
		_epsilon(getParam<Real>("epsilon")),
		_material_data(declareProperty<CLawBoundaryMaterialData>("bnd_material_data"))
{
}

void CLawBoundaryMaterial::computeProperties()
{
		_claw_problem.computeBoundaryMaterial(*this);
}

Real CLawBoundaryMaterial::penalty()
{
	Real h_face = (_current_elem_volume+_current_elem_volume)/_current_side_volume /2.;
	return _sigma*_var_order*_var_order/h_face;
}

