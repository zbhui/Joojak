
#include "CLawFaceMaterial.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawFaceMaterial>()
{
  InputParameters params = validParams<CLawMaterial>();
  params.addParam<Real>("ds", 1.490116119384766e-08, "微扰量");
  params.addParam<Real>("sigma", 6, "通量罚值，默认值为6");
  params.addParam<Real>("epsilon", 1, "对称项罚值，可以取1, 0 , -1，分别对应SIP, IIP, NIP");
  return params;
}

CLawFaceMaterial::CLawFaceMaterial(const std::string & name, InputParameters parameter):
		CLawMaterial(name, parameter),
		_ds(getParam<Real>("ds")),
		_sigma(getParam<Real>("sigma")),
		_epsilon(getParam<Real>("epsilon")),
		_material_data(declareProperty<CLawFaceMaterialData>("face_material_data"))
{
	if(_bnd && _neighbor)
	{
		for (int eq = 0; eq < _n_equations; ++eq)
		{
			MooseVariable &val = _claw_problem.getVariable(_tid, _variables[eq]);
			mooseAssert(val.order() == _var_order, "变量的阶不同");

			_ul.push_back(_is_implicit ? &val.sln() : &val.slnOld());
			_ur.push_back(_is_implicit ? &val.slnNeighbor() : &val.slnOldNeighbor());
			_grad_ul.push_back(_is_implicit ? &val.gradSln(): &val.gradSlnOld());
			_grad_ur.push_back(_is_implicit ? &val.gradSlnNeighbor(): &val.gradSlnOldNeighbor());
		}
	}
}

void CLawFaceMaterial::computeProperties()
{
	if(_bnd && _neighbor)
		_claw_problem.computeFaceMaterial(*this);
}

Real CLawFaceMaterial::penalty()
{
	Real h_face = (_current_elem_volume+_neighbor_elem_volume)/_current_side_volume /2.;
	return _sigma*_var_order*_var_order/h_face;
}
