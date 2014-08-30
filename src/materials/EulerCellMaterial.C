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

#include "EulerCellMaterial.h"

template<>
InputParameters validParams<EulerCellMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<EulerBase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

EulerCellMaterial::EulerCellMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		EulerBase(name, parameters),
		_invis_term(declareProperty<std::vector<RealVectorValue> >("cell_material")),
		_jacobi(declareProperty<std::vector<std::vector<RealVectorValue> > >("cell_jacobi_variable"))
{
	_n_equations = coupledComponents("variables");
	for (int eq = 0; eq < _n_equations; ++eq)
		_uh.push_back(&coupledValue("variables", eq));
}

void EulerCellMaterial::computeQpProperties()
{
	if(_bnd) return;

	Real uh[10];
	computeQpValue(uh);
	_invis_term[_qp].resize(_n_equations);
	inviscousTerm(_invis_term[_qp], uh);

	RealVectorValue invis_term_new[10];
	_jacobi[_qp].resize(_n_equations);
	for (int p = 0; p < _n_equations; ++p)
		_jacobi[_qp][p].resize(_n_equations);

	for (int q = 0; q < _n_equations; ++q)
	{
		uh[q] += _ds;
		inviscousTerm(invis_term_new, uh);
		for (int p = 0; p < _n_equations; ++p)
			_jacobi[_qp][p][q] = (invis_term_new[p] - _invis_term[_qp][p])/_ds;

		uh[q] -= _ds;
	}
}

void EulerCellMaterial::computeQpValue(Real* uh)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		uh[eq] = (*_uh[eq])[_qp];
}
