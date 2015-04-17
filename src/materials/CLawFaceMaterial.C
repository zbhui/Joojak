
#include "CLawFaceMaterial.h"

template<>
InputParameters validParams<CLawFaceMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addParam<Real>("ds", 1.490116119384766e-08, "微扰量");
  params.addParam<Real>("sigma", 6, "通量罚值，默认值为6");
  params.addParam<Real>("epsilon", 1, "对称项罚值，可以取1, 0 , -1，分别对应SIP, IIP, NIP");
  return params;
}

CLawFaceMaterial::CLawFaceMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		CLawInterface(parameters),
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
			MooseVariable &val = getVariable(eq);
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
		RealVectorValue lift_new[10], penalty_neighbor_new[10];
		RealVectorValue vis_term_left[10], vis_term_right[10], vis_term_new[10];

		computeQpLeftValue(ul);
		computeQpRightValue(ur);
		computeQpLeftGradValue(dul);
		computeQpRightGradValue(dur);

		liftOperator(&_lift[_qp][0], ul, ur);
		fluxTerm(&_flux[_qp][0], ul, ur, dul, dur);

		for (int q = 0; q < _n_equations; ++q)
		{
			ul[q] += _ds;
			liftOperator(lift_new, ul, ur);
			fluxTerm(flux_new, ul, ur, dul, dur);
			for (int p = 0; p < _n_equations; ++p)
			{
				Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
				_flux_jacobi_variable_ee[_qp][p][q] = tmp;
				_lift_jacobi_variable[_qp][p][q] = (lift_new[p] - _lift[_qp][p])/_ds;
			}
			ul[q] -= _ds;

			ur[q] += _ds;
			liftOperator(lift_new, ul, ur);
			fluxTerm(flux_new, ul, ur, dul, dur);
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
				fluxTerm(flux_new, ul, ur, dul, dur);
				for (int p = 0; p < _n_equations; ++p)
				{
					Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
					_flux_jacobi_grad_variable_ee[_qp][p][q](beta) = tmp;
				}
				dul[q](beta) -= _ds;

				dur[q](beta) += _ds;
				fluxTerm(flux_new, ul, ur, dul, dur);
				for (int p = 0; p < _n_equations; ++p)
				{
					Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
					_flux_jacobi_grad_variable_en[_qp][p][q](beta) = tmp;
				}
				dur[q](beta) -= _ds;
			}
		}

		addPenalty();
	}

}

void CLawFaceMaterial::computeQpLeftValue(Real* ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_ul[eq])[_qp];
}

void CLawFaceMaterial::computeQpRightValue(Real* ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_ur[eq])[_qp];
}

void CLawFaceMaterial::computeQpLeftGradValue(RealGradient *ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_grad_ul[eq])[_qp];
}

void CLawFaceMaterial::computeQpRightGradValue(RealGradient *ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_grad_ur[eq])[_qp];
}

void CLawFaceMaterial::fluxTerm(Real *flux, Real* ul, Real* ur, RealGradient *dul, RealGradient *dur)
{
	RealVectorValue ifl[5], ifr[5], vfl[5], vfr[5];

	diffusionTerm(vfl, ul, dul);
	diffusionTerm(vfr, ur, dur);

	fluxRiemann(flux, ul, ur, _normals[_qp]);
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux[eq] -= 0.5*((vfl[eq]+vfr[eq]))*_normals[_qp];
	}
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

void CLawFaceMaterial::liftOperator(RealVectorValue* lift, Real* ul, Real* ur)
{
	RealGradient duh[10];
	Real uh[10];
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		duh[eq] = (ul[eq]-ur[eq])/2.*_normals[_qp];
		uh[eq] = (ul[eq]+ur[eq])/2;
	}

	diffusionTerm(lift, uh, duh);
}

void CLawFaceMaterial::fillQpValue(QpValue &qp_value)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		qp_value.ul[eq] = (*_ul[eq])[_qp];
		qp_value.ur[eq] = (*_ur[eq])[_qp];
		qp_value.dul[eq] = (*_grad_ul[eq])[_qp];
		qp_value.dur[eq] = (*_grad_ur[eq])[_qp];
		qp_value.duh[eq] = (qp_value.ul[eq]-qp_value.ur[eq])/2.*_normals[_qp];
		qp_value.u_bar[eq] = (qp_value.ul[eq]+qp_value.ur[eq])/2;
	}
}

void CLawFaceMaterial::computeQpLift(QpValue &qp_value)
{
	diffusionTerm(qp_value.lift, qp_value.u_bar, qp_value.duh);
}

void CLawFaceMaterial::computeQpLift(RealVectorValue* lift, QpValue& qp_value)
{
	diffusionTerm(lift, qp_value.u_bar, qp_value.duh);
}

void CLawFaceMaterial::computeQpFlux(Real* flux, QpValue& qp_value)
{
	RealVectorValue ifl[5], ifr[5], vfl[5], vfr[5];

	diffusionTerm(vfl, qp_value.ul, qp_value.dul);
	diffusionTerm(vfr, qp_value.ur, qp_value.dur);

	fluxRiemann(flux, qp_value.ul, qp_value.ur, _normals[_qp]);
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux[eq] -= 0.5*((vfl[eq]+vfr[eq]))*_normals[_qp];
	}
}

void CLawFaceMaterial::addPenalty()
{

	Real h_face = (_current_elem_volume+_neighbor_elem_volume)/_current_side_volume /2.;
	Real penalty = _sigma*_var_order*_var_order/h_face;
	for (int p = 0; p < _n_equations; ++p)
	{
		_flux[_qp][p] += penalty*_lift[_qp][p]*_normals[_qp];
		for (int q = 0; q < _n_equations; ++q)
		{
			_flux_jacobi_variable_ee[_qp][p][q] += penalty*_lift_jacobi_variable[_qp][p][q]*_normals[_qp];
			_flux_jacobi_variable_en[_qp][p][q] += penalty*_lift_jacobi_variable_neighbor[_qp][p][q]*_normals[_qp];
		}
	}
}
