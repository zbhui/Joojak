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
  params += validParams<CFDBase>();

  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

EulerCellMaterial::EulerCellMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		CFDBase(name, parameters),
		_invis_term(declareProperty<std::vector<RealVectorValue> >("invis_term_mt"))
{
	_n_equations = coupledComponents("variables");
	for (int eq = 0; eq < _n_equations; ++eq)
		_uh.push_back(&coupledValue("varialbes", eq));

	for (int qp = 0; qp < _qrule->n_points(); ++qp)
	{
		_invis_term[qp].resize(_n_equations);
	}
}

void EulerCellMaterial::computeQpProperties()
{
	Real uh[10];
	computeQpValue(uh);
	inviscousTerm(_invis_term[_qp], uh);
}

void EulerCellMaterial::computeQpValue(Real* uh)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		uh[eq] = (*_uh[eq])[_qp];
}
