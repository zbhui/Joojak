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
#include "CLawProblem.h"

template<>
InputParameters validParams<IsoVortexBndMaterial>()
{
  InputParameters params = validParams<CLawBoundaryMaterial>();
  params += validParams<IsoVortexBase>();
  params.set<std::string>("bc_type") = "none";
  return params;
}

IsoVortexBndMaterial::IsoVortexBndMaterial(const std::string & name, InputParameters parameters):
	CLawBoundaryMaterial(name, parameters),
	IsoVortexBase(name, parameters)
{
}

void IsoVortexBndMaterial::computeQpRightValue(Real* ur)
{
	Point p = _q_point[_qp];
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = IsoVortexBase::value(_t, p, eq);
}

void IsoVortexBndMaterial::computeQpFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul)
{
	Real h_face = (_current_elem_volume+_current_elem_volume)/_current_side_volume /2.;
	Real penalty = _sigma*_var_order*_var_order/h_face;
	Point normal = _normals[_qp];

	Real ur[10];
	RealGradient dur[10];

	for (int eq = 0; eq < _n_equations; ++eq)
		dur[eq] = dul[eq];

	computeQpRightValue(ur);
	_claw_problem.computeFaceFlux(flux, lift, ul, ur, dul, dur, normal, penalty);
}
