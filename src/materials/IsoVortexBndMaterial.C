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

#include "IsoVortexBndMaterial.h"

template<>
InputParameters validParams<IsoVortexBndMaterial>()
{
  InputParameters params = validParams<EulerBndMaterial>();
  params += validParams<IsoVortexBase>();
  return params;
}

IsoVortexBndMaterial::IsoVortexBndMaterial(const std::string & name, InputParameters parameters):
		EulerBndMaterial(name, parameters),
		IsoVortexBase(name, parameters)
{
}

void IsoVortexBndMaterial::computeQpRightValue(Real* ur)
{
	Point p = _q_point[_qp];
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = IsoVortexBase::value(_t, p, eq);
}
