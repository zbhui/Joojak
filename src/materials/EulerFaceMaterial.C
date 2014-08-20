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
		_invis_term(declareProperty<std::vector<RealVectorValue> >("left_material")),
		_invis_term_neighbor(declareProperty<std::vector<RealVectorValue> >("right_material")),
		_flux_diff(declareProperty<Real>("flux_diff"))
{
	_n_equations = coupledComponents("variables");
	for (size_t eq = 0; eq < _n_equations; ++eq)
	{
		MooseVariable &val = *getVar("variables", eq);
		_ul.push_back(_is_implicit ? &val.sln() : &val.slnOld());
		_ur.push_back(_is_implicit ? &val.slnNeighbor() : &val.slnOldNeighbor());
	}
}

void EulerFaceMaterial::computeQpProperties()
{
	Real ul[10], ur[10];
	computeQpLeftValue(ul);
	computeQpRightValue(ur);

	_invis_term[_qp].resize(_n_equations);
	_invis_term_neighbor[_qp].resize(_n_equations);

	inviscousTerm(_invis_term[_qp], ul);
	inviscousTerm(_invis_term_neighbor[_qp], ur);
	_flux_diff[_qp] = 1.;
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
