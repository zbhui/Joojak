
#include "CLawBoundaryMaterial.h"

template<>
InputParameters validParams<CLawBoundaryMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addParam<Real>("ds", 1.490116119384766e-08, "微扰量");
  params.addParam<Real>("sigma", 6, "通量罚值，默认值为6");
  params.addParam<Real>("epsilon", 1, "对称项罚值，可以取1, 0 , -1，分别对应SIP, IIP, NIP");
  return params;
}

CLawBoundaryMaterial::CLawBoundaryMaterial(const std::string & name, InputParameters parameters):
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
		_lift_jacobi_variable_neighbor(declareProperty<std::vector<std::vector<RealVectorValue> > >("lift_jacobi_variable_en"))
{
	if(_bnd)
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

void CLawBoundaryMaterial::computeQpProperties()
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

void CLawBoundaryMaterial::computeQpLeftValue(Real* ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_ul[eq])[_qp];
}

void CLawBoundaryMaterial::computeQpRightValue(Real* ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_ur[eq])[_qp];
}

void CLawBoundaryMaterial::computeQpLeftGradValue(RealGradient *ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_grad_ul[eq])[_qp];
}

void CLawBoundaryMaterial::computeQpRightGradValue(RealGradient *ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_grad_ul[eq])[_qp];
}

void CLawBoundaryMaterial::fluxTerm(Real *flux, Real* ul, Real* ur, RealGradient *dul, RealGradient *dur)
{
	RealVectorValue ifl[5], ifr[5], vfl[5], vfr[5];

	convertionTerm(ifl, ul);
	convertionTerm(ifr, ur);
	diffusionTerm(vfl, ul, dul);
	diffusionTerm(vfr, ur, dur);

	Real lam = 1;//(maxEigenValue(ul, _normals[_qp]) + maxEigenValue(ur, _normals[_qp]))/2.;
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux[eq] = 0.5*(ifl[eq] + ifr[eq] - (vfl[eq]+vfr[eq]))*_normals[_qp] + lam*(ul[eq] - ur[eq]);
	}
}

void CLawBoundaryMaterial::resizeQpProperty()
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

void CLawBoundaryMaterial::liftOperator(RealVectorValue* lift, Real* ul, Real* ur)
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

void CLawBoundaryMaterial::computeQpValue(Real* ul, Real* ur, RealGradient* dul, RealGradient* dur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		ul[eq] = (*_ul[eq])[_qp];
		dul[eq] = (*_grad_ul[eq])[_qp];

		dur[eq] = (*_grad_ul[eq])[_qp];
	}
}

void CLawBoundaryMaterial::addPenalty()
{

	Real h_face = (_current_elem_volume+_current_elem_volume)/_current_side_volume /2.;
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
