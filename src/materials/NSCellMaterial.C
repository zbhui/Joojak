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

#include "NSCellMaterial.h"

template<>
InputParameters validParams<NSCellMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<NSBase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

NSCellMaterial::NSCellMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		NSBase(name, parameters),
		_flux_term(declareProperty<std::vector<RealVectorValue> >("flux_term")),
		_jacobi_variable(declareProperty<std::vector<std::vector<RealVectorValue> > >("jacobi_variable")),
		_jacobi_grad_variable(declareProperty<std::vector<std::vector<RealTensorValue> > >("jacobi_grad_variable"))
{
	_n_equations = coupledComponents("variables");
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		_uh.push_back(&coupledValue("variables", eq));
		_grad_uh.push_back(&coupledGradient("variables", eq));
	}
}

void NSCellMaterial::computeQpProperties()
{
	if(_bnd) return;

	resizeQpProperty();

	Real uh[10];
	RealGradient duh[10];
	computeQpValue(uh, duh);
	fluxTerm(&_flux_term[_qp][0], uh, duh);

//	RealVectorValue invis_term_new[10];
//	for (int q = 0; q < _n_equations; ++q)
//	{
//		uh[q] += _ds;
//		inviscousTerm(invis_term_new, uh);
//		for (int p = 0; p < _n_equations; ++p)
//			_jacobi_variable[_qp][p][q] = (invis_term_new[p] - _invis_term[_qp][p])/_ds;
//
//		uh[q] -= _ds;
//	}
}

void NSCellMaterial::computeQpValue(Real* uh, RealGradient *duh)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		uh[eq] = (*_uh[eq])[_qp];
		duh[eq] = (*_grad_uh[eq])[_qp];
	}
}

void NSCellMaterial::resizeQpProperty()
{
	_flux_term[_qp].resize(_n_equations);
	_jacobi_variable[_qp].resize(_n_equations);
	_jacobi_grad_variable[_qp].resize(_n_equations);
	for (int p = 0; p < _n_equations; ++p)
	{
		_jacobi_variable[_qp][p].resize(_n_equations);
		_jacobi_grad_variable[_qp][p].resize(_n_equations);

	}

}

void NSCellMaterial::fluxTerm(RealVectorValue* flux_term, Real* uh, RealGradient* duh)
{
	RealVectorValue invis_term[10], viscous_term[10];
	inviscousTerm(invis_term, uh);
	viscousTerm(viscous_term, uh, duh);
	for (int eq = 0; eq < _n_equations; ++eq)
		flux_term[eq] = invis_term[eq] - viscous_term[eq];
}
