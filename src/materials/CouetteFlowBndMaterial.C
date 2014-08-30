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

#include "CouetteFlowBndMaterial.h"

template<>
InputParameters validParams<CouetteFlowBndMaterial>()
{
  InputParameters params = validParams<NSBndMaterial>();
  params += validParams<CouetteFlowBase>();
  return params;
}

CouetteFlowBndMaterial::CouetteFlowBndMaterial(const std::string & name, InputParameters parameters):
		NSBndMaterial(name, parameters),
		CouetteFlowBase(name, parameters)
{
}

void CouetteFlowBndMaterial::computeQpRightValue(Real *ur,RealGradient *dur, Real *ul, RealGradient *dul)
{
	Point p = _q_point[_qp];
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		ur[eq] = CouetteFlowBase::value(_t, p, eq);
		dur[eq] = dul[eq];
	}
}
