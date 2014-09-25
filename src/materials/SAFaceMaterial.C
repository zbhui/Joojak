/** ***********************************************************
*  @file
*  @brief
*  @author	刘  伟
*
*  Program:   Joojak
*  Copyright (c) 刘伟，张来平，2014，空气动力学国家重点实验室(SKLA)
*  All rights reserved.
*  ************************************************************
**/

#include "SAFaceMaterial.h"

template<>
InputParameters validParams<SAFaceMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<SABase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

SAFaceMaterial::SAFaceMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		SABase(name, parameters),
		_flux(declareProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable_ee(declareProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_ee")),
		_flux_jacobi_variable_en(declareProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_en")),
		_flux_jacobi_variable_ne(declareProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_ne")),
		_flux_jacobi_variable_nn(declareProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_nn")),
		_flux_jacobi_grad_variable_ee(declareProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_ee")),
		_flux_jacobi_grad_variable_en(declareProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_en")),
		_flux_jacobi_grad_variable_ne(declareProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_ne")),
		_flux_jacobi_grad_variable_nn(declareProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_nn")),

		_penalty(declareProperty<std::vector<RealVectorValue> >("penalty")),
		_penalty_neighbor(declareProperty<std::vector<RealVectorValue> >("penalty_neighbor")),
		_penalty_jacobi_variable_ee(declareProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_ee")),
		_penalty_jacobi_variable_en(declareProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_en")),
		_penalty_jacobi_variable_ne(declareProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_ne")),
		_penalty_jacobi_variable_nn(declareProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_nn"))

{
	_n_equations = coupledComponents("variables");

	if(_bnd && _neighbor)
	{
		for (int eq = 0; eq < _n_equations; ++eq)
		{
			MooseVariable &val = *getVar("variables", eq);
			_ul.push_back(_is_implicit ? &val.sln() : &val.slnOld());
			_ur.push_back(_is_implicit ? &val.slnNeighbor() : &val.slnOldNeighbor());
			_grad_ul.push_back(_is_implicit ? &val.gradSln(): &val.gradSlnOld());
			_grad_ur.push_back(_is_implicit ? &val.gradSlnNeighbor(): &val.gradSlnOldNeighbor());
		}
	}
}

void SAFaceMaterial::computeQpProperties()
{
	if(_bnd && _neighbor)
	{
		resizeQpProperty();

		Real ul[10], ur[10], ur_new[10];
		RealGradient dul[10], dur[10];
		Real flux_new[10];
		RealVectorValue penalty_new[10], penalty_neighbor_new[10];
		RealVectorValue vis_term_left[10], vis_term_right[10], vis_term_new[10];

		computeQpLeftValue(ul);
		computeQpRightValue(ur);
		computeQpLeftGradValue(dul);
		computeQpRightGradValue(dur);

		penaltyTerm(&_penalty[_qp][0], &_penalty_neighbor[_qp][0], ul, ur);
		fluxTerm(&_flux[_qp][0], ul, ur, dul, dur);

		for (int q = 0; q < _n_equations; ++q)
		{
			ul[q] += _ds;
			penaltyTerm(penalty_new, penalty_neighbor_new, ul, ur);
			fluxTerm(flux_new, ul, ur, dul, dur);
			for (int p = 0; p < _n_equations; ++p)
			{
				Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
				_flux_jacobi_variable_ee[_qp][p][q] = tmp;
				_flux_jacobi_variable_ne[_qp][p][q] = -tmp;

				_penalty_jacobi_variable_ee[_qp][p][q] = (penalty_new[p] - _penalty[_qp][p])/_ds;
				_penalty_jacobi_variable_ne[_qp][p][q] = (penalty_neighbor_new[p] - _penalty_neighbor[_qp][p])/_ds;
			}
			ul[q] -= _ds;

			ur[q] += _ds;
			penaltyTerm(penalty_new, penalty_neighbor_new, ul, ur);
			fluxTerm(flux_new, ul, ur, dul, dur);
			for (int p = 0; p < _n_equations; ++p)
			{
				Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
				_flux_jacobi_variable_en[_qp][p][q] = tmp;
				_flux_jacobi_variable_nn[_qp][p][q] = -tmp;

				_penalty_jacobi_variable_en[_qp][p][q] = (penalty_new[p] - _penalty[_qp][p])/_ds;
				_penalty_jacobi_variable_nn[_qp][p][q] = (penalty_neighbor_new[p] - _penalty_neighbor[_qp][p])/_ds;
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
					_flux_jacobi_grad_variable_ne[_qp][p][q](beta) = -tmp;
				}
				dul[q](beta) -= _ds;

				dur[q](beta) += _ds;
				fluxTerm(flux_new, ul, ur, dul, dur);
				for (int p = 0; p < _n_equations; ++p)
				{
					Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
					_flux_jacobi_grad_variable_en[_qp][p][q](beta) = tmp;
					_flux_jacobi_grad_variable_nn[_qp][p][q](beta) = -tmp;
				}
				dur[q](beta) -= _ds;
			}
		}

	}
}

void SAFaceMaterial::computeQpLeftValue(Real* ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_ul[eq])[_qp];
}

void SAFaceMaterial::computeQpRightValue(Real* ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_ur[eq])[_qp];
}

void SAFaceMaterial::computeQpLeftGradValue(RealGradient *ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_grad_ul[eq])[_qp];
}

void SAFaceMaterial::computeQpRightGradValue(RealGradient *ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_grad_ur[eq])[_qp];
}

void SAFaceMaterial::fluxTerm(Real *flux, Real* ul, Real* ur, RealGradient *dul, RealGradient *dur)
{
	RealVectorValue ifl[10], ifr[10], vfl[10], vfr[10];

	inviscousTerm(ifl, ul);
	inviscousTerm(ifr, ur);
	viscousTerm(vfl, ul, dul);
	viscousTerm(vfr, ur, dur);

	Real lam = (maxEigenValue(ul, _normals[_qp]) + maxEigenValue(ur, _normals[_qp]))/2.;
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux[eq] = 0.5*(ifl[eq] + ifr[eq] - (vfl[eq]+vfr[eq]))*_normals[_qp] + lam*(ul[eq] - ur[eq]);
	}
}

void SAFaceMaterial::penaltyTerm(RealVectorValue* penalty, RealVectorValue* penalty_neighbor, Real* ul, Real* ur)
{
	RealGradient duh[10];
	for (int eq = 0; eq < _n_equations; ++eq)
		duh[eq] = (ul[eq]-ur[eq])/2.*_normals[_qp];

	viscousTerm(penalty, ul, duh);
	viscousTerm(penalty_neighbor, ur, duh);
}

void SAFaceMaterial::resizeQpProperty()
{
	_flux[_qp].resize(_n_equations);
	_flux_jacobi_variable_ee[_qp].resize(_n_equations);
	_flux_jacobi_variable_en[_qp].resize(_n_equations);
	_flux_jacobi_variable_ne[_qp].resize(_n_equations);
	_flux_jacobi_variable_nn[_qp].resize(_n_equations);
	_flux_jacobi_grad_variable_ee[_qp].resize(_n_equations);
	_flux_jacobi_grad_variable_en[_qp].resize(_n_equations);
	_flux_jacobi_grad_variable_ne[_qp].resize(_n_equations);
	_flux_jacobi_grad_variable_nn[_qp].resize(_n_equations);

	_penalty[_qp].resize(_n_equations);
	_penalty_neighbor[_qp].resize(_n_equations);
	_penalty_jacobi_variable_ee[_qp].resize(_n_equations);
	_penalty_jacobi_variable_en[_qp].resize(_n_equations);
	_penalty_jacobi_variable_ne[_qp].resize(_n_equations);
	_penalty_jacobi_variable_nn[_qp].resize(_n_equations);

	for (int p = 0; p < _n_equations; ++p)
	{
		_flux_jacobi_variable_ee[_qp][p].resize(_n_equations);
		_flux_jacobi_variable_en[_qp][p].resize(_n_equations);
		_flux_jacobi_variable_ne[_qp][p].resize(_n_equations);
		_flux_jacobi_variable_nn[_qp][p].resize(_n_equations);
		_flux_jacobi_grad_variable_ee[_qp][p].resize(_n_equations);
		_flux_jacobi_grad_variable_en[_qp][p].resize(_n_equations);
		_flux_jacobi_grad_variable_ne[_qp][p].resize(_n_equations);
		_flux_jacobi_grad_variable_nn[_qp][p].resize(_n_equations);

		_penalty_jacobi_variable_ee[_qp][p].resize(_n_equations);
		_penalty_jacobi_variable_en[_qp][p].resize(_n_equations);
		_penalty_jacobi_variable_ne[_qp][p].resize(_n_equations);
		_penalty_jacobi_variable_nn[_qp][p].resize(_n_equations);
	}
}


