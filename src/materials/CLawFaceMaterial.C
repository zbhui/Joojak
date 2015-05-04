
#include "CLawFaceMaterial.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawFaceMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addParam<Real>("ds", 1.490116119384766e-08, "微扰量");
  params.addParam<Real>("sigma", 6, "通量罚值，默认值为6");
  params.addParam<Real>("epsilon", 1, "对称项罚值，可以取1, 0 , -1，分别对应SIP, IIP, NIP");
  return params;
}

CLawFaceMaterial::CLawFaceMaterial(const std::string & name, InputParameters parameter):
		Material(name, parameter),
		_claw_problem(static_cast<CLawProblem&>(_fe_problem)),
		_nl(_claw_problem.getNonlinearSystem()),
		_tid(parameter.get<THREAD_ID>("_tid")),
		_variables(_nl.getVariableNames()),
		_n_equations(_variables.size()),
		_var_order(_claw_problem.getVariable(_tid, _variables[0]).order()),

		_current_elem_volume(_assembly.elemVolume()),
		_neighbor_elem_volume(_assembly.neighborVolume()),
		_current_side_volume(_assembly.sideElemVolume()),
		_ds(getParam<Real>("ds")),
		_sigma(getParam<Real>("sigma")),
		_epsilon(getParam<Real>("epsilon")),
		_flux(declareProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable_ee(declareProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_ee")),
		_flux_jacobi_variable_en(declareProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_en")),
		_flux_jacobi_grad_variable_ee(declareProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_ee")),
		_flux_jacobi_grad_variable_en(declareProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_en")),

		_lift(declareProperty<std::vector<RealVectorValue> >("lift")),
		_lift_jacobi_variable(declareProperty<std::vector<std::vector<RealVectorValue> > >("lift_jacobi_variable_ee")),
		_lift_jacobi_variable_neighbor(declareProperty<std::vector<std::vector<RealVectorValue> > >("lift_jacobi_variable_en")),
		_face_material_data(declareProperty<CLawFaceMaterialData>("face_material_data"))
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

void CLawFaceMaterial::computeQpProperties()
{
	if(_bnd && _neighbor)
	{
		resizeQpProperty();

		Real ul[10], ur[10], ur_new[10];
		RealGradient dul[10], dur[10];
		Real flux_new[10];
		RealVectorValue lift_new[10];
		RealVectorValue vis_term_left[10], vis_term_right[10], vis_term_new[10];

		computeQpValue(ul, ur, dul, dur);

		computeQpFlux(&_flux[_qp][0], &_lift[_qp][0], ul, ur, dul, dur);
		for (int q = 0; q < _n_equations; ++q)
		{
			ul[q] += _ds;
			computeQpFlux(flux_new, lift_new, ul, ur, dul, dur);
			for (int p = 0; p < _n_equations; ++p)
			{
				Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
				_flux_jacobi_variable_ee[_qp][p][q] = tmp;
				_lift_jacobi_variable[_qp][p][q] = (lift_new[p] - _lift[_qp][p])/_ds;
			}
			ul[q] -= _ds;

			ur[q] += _ds;
			computeQpFlux(flux_new, lift_new, ul, ur, dul, dur);
			for (int p = 0; p < _n_equations; ++p)
			{
				Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
				_flux_jacobi_variable_en[_qp][p][q] = tmp;
				_lift_jacobi_variable_neighbor[_qp][p][q] = (lift_new[p] - _lift[_qp][p])/_ds;
			}
			ur[q] -= _ds;

			for (int beta = 0; beta < 3; ++beta)
			for (int q = 0; q < _n_equations; ++q)
			{
				dul[q](beta) += _ds;
				computeQpFlux(flux_new, lift_new, ul, ur, dul, dur);
				for (int p = 0; p < _n_equations; ++p)
				{
					Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
					_flux_jacobi_grad_variable_ee[_qp][p][q](beta) = tmp;
				}
				dul[q](beta) -= _ds;

				dur[q](beta) += _ds;
				computeQpFlux(flux_new, lift_new, ul, ur, dul, dur);
				for (int p = 0; p < _n_equations; ++p)
				{
					Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
					_flux_jacobi_grad_variable_en[_qp][p][q](beta) = tmp;
				}
				dur[q](beta) -= _ds;
			}
		}
	}

}

void CLawFaceMaterial::computeQpValue(Real *ul, Real *ur, RealGradient *dul, RealGradient *dur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		ul[eq] = (*_ul[eq])[_qp];
		ur[eq] = (*_ur[eq])[_qp];
		dul[eq] = (*_grad_ul[eq])[_qp];
		dur[eq] = (*_grad_ur[eq])[_qp];
	}
}

void CLawFaceMaterial::computeQpFlux(Real* flux, RealVectorValue* lift, Real* ul, Real* ur, RealGradient* dul, RealGradient* dur)
{
	Real h_face = (_current_elem_volume+_neighbor_elem_volume)/_current_side_volume /2.;
	Real penalty = _sigma*_var_order*_var_order/h_face;
	Point normal = _normals[_qp];
	_claw_problem.computeFaceFlux(flux, lift, ul, ur, dul, dur, normal, penalty);
}

void CLawFaceMaterial::resizeQpProperty()
{
	_flux[_qp].resize(_n_equations);
	_flux_jacobi_variable_ee[_qp].resize(_n_equations);
	_flux_jacobi_variable_en[_qp].resize(_n_equations);
	_flux_jacobi_grad_variable_ee[_qp].resize(_n_equations);
	_flux_jacobi_grad_variable_en[_qp].resize(_n_equations);

	_lift[_qp].resize(_n_equations);
	_lift_jacobi_variable[_qp].resize(_n_equations);
	_lift_jacobi_variable_neighbor[_qp].resize(_n_equations);

	for (int p = 0; p < _n_equations; ++p)
	{
		_flux_jacobi_variable_ee[_qp][p].resize(_n_equations);
		_flux_jacobi_variable_en[_qp][p].resize(_n_equations);
		_flux_jacobi_grad_variable_ee[_qp][p].resize(_n_equations);
		_flux_jacobi_grad_variable_en[_qp][p].resize(_n_equations);

		_lift_jacobi_variable[_qp][p].resize(_n_equations);
		_lift_jacobi_variable_neighbor[_qp][p].resize(_n_equations);
	}
}

