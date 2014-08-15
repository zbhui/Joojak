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

#include "EulerMaterial.h"

template<>
InputParameters validParams<EulerMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  params.addRequiredParam<Real>("mach",     "马赫数");
  params.addRequiredParam<Real>("reynolds", "雷诺数");
  params.addParam<Real>("gamma", 1.4, "比热比");
  params.addParam<Real>("prandtl", 0.72, "prandtl数");

  return params;
}

EulerMaterial::EulerMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
	    _gamma(getParam<Real>("gamma")),
	    _prandtl(getParam<Real>("prandtl")),
	    _reynolds(getParam<Real>("reynolds")),
	    _mach(getParam<Real>("mach")),
	    _inviscous_term(declareProperty<RealVectorValue*>("inviscous"))
{
	_n_equations = coupledComponents("variables");
	for (int eq = 0; eq < _n_equations; ++eq)
		_uh.push_back(&coupledValue("varialbes", eq));
}

void EulerMaterial::computeQpProperties()
{
	Real uh[10];
	computeQpValue(uh);

	Real rho, p, h;
	Real u, v, w;
	rho = uh[0];
	u = uh[1]/rho;
	v = uh[2]/rho;
	w = uh[3]/rho;
	p = pressure(uh);
	h = enthalpy(uh);

	int eq = 0;

	eq = 0;
	_inviscous_term[_qp][eq](0) = uh[1];	// rhou
	_inviscous_term[_qp][eq](1) = uh[2];	// rhov
	_inviscous_term[_qp][eq](2) = uh[3];	// rhow

	eq = 1;
	_inviscous_term[_qp][eq](0) = uh[1] * u + p;
	_inviscous_term[_qp][eq](1) = uh[1] * v;
	_inviscous_term[_qp][eq](2) = uh[1] * w;

	eq = 2;
	_inviscous_term[_qp][eq](0) = uh[2] * u;
	_inviscous_term[_qp][eq](1) = uh[2] * v + p;
	_inviscous_term[_qp][eq](2) = uh[2] * w;

	eq = 3;
	_inviscous_term[_qp][eq](0) = uh[3] * u;
	_inviscous_term[_qp][eq](1) = uh[3] * v;
	_inviscous_term[_qp][eq](2) = uh[3] * w + p;

	eq = 4;
	_inviscous_term[_qp][eq](0) = rho * h * u;
	_inviscous_term[_qp][eq](1) = rho * h * v;
	_inviscous_term[_qp][eq](2) = rho * h * w;
}

void EulerMaterial::computeQpValue(Real* uh)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		uh[eq] = (*_uh[eq])[_qp];
}

Real EulerMaterial::pressure(Real *uh)
{
	return (_gamma-1) * (uh[4] - 0.5*(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0]);  //
}
Real EulerMaterial::physicalViscosity(Real* uh)
{
	return 1.0;
}

Real EulerMaterial::enthalpy(Real *uh)
{
	Real p = pressure(uh);
	return (uh[4] + p)/uh[0];
}

Real EulerMaterial::temperature(Real* uh)
{
	Real p = pressure(uh);
	return _gamma*_mach*_mach*p/uh[0];
}

Real EulerMaterial::mach_local(Real* uh)
{
	Real vel = std::sqrt(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0];
	Real c = std::sqrt(temperature(uh))/_mach;
	return vel/c;
}
