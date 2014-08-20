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
  params += validParams<CFDBase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

EulerFaceMaterial::EulerFaceMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		CFDBase(name, parameters),
		_flux(declareProperty<std::vector<Real> >("flux"))
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
		RealVectorValue invis_term[10], invis_term_neighbor[10];

		computeQpLeftValue(ul);
		computeQpRightValue(ur);
		inviscousTerm(invis_term, ul);
		inviscousTerm(invis_term_neighbor, ur);

		Real lam = (maxEigenValue(ul, _normals[_qp]) + maxEigenValue(ur, _normals[_qp]))/2.;
		for (int eq = 0; eq < _n_equations; ++eq)
			_flux[_qp][eq] = 0.5*(invis_term[eq] + invis_term_neighbor[eq])*_normals[_qp] + lam*(ul[eq]-ur[eq]);
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
