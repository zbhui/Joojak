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

#include "EulerFaceMaterial.h"

template<>
InputParameters validParams<EulerFaceMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<EulerBase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

EulerFaceMaterial::EulerFaceMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		EulerBase(name, parameters),
		_flux(declareProperty<std::vector<Real> >("flux")),
		_jacobi_variable_ee(declareProperty<std::vector<std::vector<Real> > >("face_jacobi_ee")),
		_jacobi_variable_en(declareProperty<std::vector<std::vector<Real> > >("face_jacobi_en")),
		_jacobi_variable_ne(declareProperty<std::vector<std::vector<Real> > >("face_jacobi_ne")),
		_jacobi_variable_nn(declareProperty<std::vector<std::vector<Real> > >("face_jacobi_nn"))
{
	_n_equations = coupledComponents("variables");

	if(_bnd && _neighbor)
	{
		for (int eq = 0; eq < _n_equations; ++eq)
		{
			MooseVariable &val = *getVar("variables", eq);
			_ul.push_back(_is_implicit ? &val.sln() : &val.slnOld());
			_ur.push_back(_is_implicit ? &val.slnNeighbor() : &val.slnOldNeighbor());
		}
	}
}

void EulerFaceMaterial::computeQpProperties()
{
	if(_bnd && _neighbor)
	{
		_flux[_qp].resize(_n_equations);

		Real ul[10], ur[10];
		Real flux_new[10];

		computeQpLeftValue(ul);
		computeQpRightValue(ur);
		fluxRiemann(&_flux[_qp][0], ul, ur);

		_jacobi_variable_ee[_qp].resize(_n_equations);
		_jacobi_variable_en[_qp].resize(_n_equations);
		_jacobi_variable_ne[_qp].resize(_n_equations);
		_jacobi_variable_nn[_qp].resize(_n_equations);
		for (int p = 0; p < _n_equations; ++p)
		{
			_jacobi_variable_ee[_qp][p].resize(_n_equations);
			_jacobi_variable_en[_qp][p].resize(_n_equations);
			_jacobi_variable_ne[_qp][p].resize(_n_equations);
			_jacobi_variable_nn[_qp][p].resize(_n_equations);
		}

		for (int q = 0; q < _n_equations; ++q)
		{
			ul[q] += _ds;
			fluxRiemann(flux_new, ul, ur);
			for (int p = 0; p < _n_equations; ++p)
			{
				Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
				_jacobi_variable_ee[_qp][p][q] = tmp;
				_jacobi_variable_ne[_qp][p][q] = -tmp;
			}
			ul[q] -= _ds;

			ur[q] += _ds;
			fluxRiemann(flux_new, ul, ur);
			for (int p = 0; p < _n_equations; ++p)
			{
				Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
				_jacobi_variable_en[_qp][p][q] = tmp;
				_jacobi_variable_nn[_qp][p][q] = -tmp;
			}
			ur[q] -= _ds;
		}
	}
}

void EulerFaceMaterial::computeQpLeftValue(Real* ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_ul[eq])[_qp];
}

void EulerFaceMaterial::computeQpRightValue(Real* ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_ur[eq])[_qp];
}

void EulerFaceMaterial::fluxRiemann(Real *flux, Real* ul, Real* ur)
{
	const Point &normal = _normals[_qp];
	RealVectorValue fl[10], fr[10];
	inviscousTerm(fl, ul);
	inviscousTerm(fr, ur);

	Real lam = (maxEigenValue(ul, _normals[_qp]) + maxEigenValue(ur, _normals[_qp]))/2.;
	for (int eq = 0; eq < _n_equations; ++eq)
		flux[eq] = 0.5*(fl[eq] + fr[eq])*normal + lam*(ul[eq] - ur[eq]);
}
