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

#include "NSFaceMaterial.h"

template<>
InputParameters validParams<NSFaceMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<NSBase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

NSFaceMaterial::NSFaceMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		NSBase(name, parameters),
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

void NSFaceMaterial::computeQpProperties()
{
	if(_bnd && _neighbor)
	{
		resizeQpProperty();

		Real ul[10], ur[10], ur_new[10];
		RealGradient dul[10], dur[10], duh[10];
		Real flux_new[10];

		computeQpLeftValue(ul);
		computeQpRightValue(ur);
		computeQpLeftGradValue(dul);
		computeQpRightGradValue(dur);

		for (int eq = 0; eq < _n_equations; ++eq)
			duh[eq] = (ul[eq]-ur[eq])*_normals[_qp];

		viscousTerm(&_penalty[_qp][0], ul, duh);
		viscousTerm(&_penalty_neighbor[_qp][0], ur, duh);

		fluxRiemann(&_flux[_qp][0], ul, ur, dul, dur);

	}
}

void NSFaceMaterial::computeQpLeftValue(Real* ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_ul[eq])[_qp];
}

void NSFaceMaterial::computeQpRightValue(Real* ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_ur[eq])[_qp];
}

void NSFaceMaterial::computeQpLeftGradValue(RealGradient *ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_grad_ul[eq])[_qp];
}

void NSFaceMaterial::computeQpRightGradValue(RealGradient *ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_grad_ur[eq])[_qp];
}

void NSFaceMaterial::fluxRiemann(Real *flux, Real* ul, Real* ur, RealGradient *dul, RealGradient *dur)
{
	RealVectorValue ifl[5], ifr[5], vfl[5], vfr[5];

	inviscousTerm(ifl, ul);
	inviscousTerm(ifr, ur);
	viscousTerm(vfl, ul, dul);
	viscousTerm(vfr, ur, dur);

	Real lam = (maxEigenValue(ul, _normals[_qp]) + maxEigenValue(ur, _normals[_qp]))/2.;
	for (int eq = 0; eq < _n_equations; ++eq)
		flux[eq] = 0.5*(ifl[eq] + ifr[eq] - (vfl[eq]+vfr[eq]))*_normals[_qp] + lam*(ul[eq] - ur[eq]);
}

void NSFaceMaterial::resizeQpProperty()
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
