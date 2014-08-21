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

#include "EulerBndMaterial.h"

template<>
InputParameters validParams<EulerBndMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<CFDBase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

EulerBndMaterial::EulerBndMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		CFDBase(name, parameters),
		_flux(declareProperty<std::vector<Real> >("flux")),
		_jacobi_variable(declareProperty<std::vector<std::vector<Real> > >("bnd_jacobi_variable"))
{
	_n_equations = coupledComponents("variables");
	for (size_t eq = 0; eq < _n_equations; ++eq)
	{
		MooseVariable &val = *getVar("variables", eq);
		_ul.push_back(_is_implicit ? &val.sln() : &val.slnOld());
	}
}

void EulerBndMaterial::computeQpProperties()
{
	_flux[_qp].resize(_n_equations);

	Real ul[10], ur[10];
	Real flux_new[10];

	computeQpLeftValue(ul);
	computeQpRightValue(ur);

	fluxRiemann(&_flux[_qp][0], ul, ur);

	_jacobi_variable[_qp].resize(_n_equations);
	for (int p = 0; p < _n_equations; ++p)
		_jacobi_variable[_qp][p].resize(_n_equations);

	for (int q = 0; q < _n_equations; ++q)
	{
		ul[q] += _ds;
		fluxRiemann(flux_new, ul, ur);
		for (int p = 0; p < _n_equations; ++p)
			_jacobi_variable[_qp][p][q] = (flux_new[p] - _flux[_qp][p])/_ds;

		ul[q] -= _ds;
	}
}

void EulerBndMaterial::computeQpLeftValue(Real* ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_ul[eq])[_qp];
}

void EulerBndMaterial::computeQpRightValue(Real* ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_ul[eq])[_qp];
}

void EulerBndMaterial::fluxRiemann(Real *flux, Real* ul, Real* ur)
{
	const Point &normal = _normals[_qp];
	RealVectorValue fl[10], fr[10];
	inviscousTerm(fl, ul);
	inviscousTerm(fr, ur);

	Real lam = (maxEigenValue(ul, _normals[_qp]) + maxEigenValue(ur, _normals[_qp]))/2.;
	for (int eq = 0; eq < _n_equations; ++eq)
		flux[eq] = 0.5*(fl[eq] + fr[eq])*normal + lam*(ul[eq] - ur[eq]);
}
