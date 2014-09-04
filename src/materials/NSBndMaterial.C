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

#include "NSBndMaterial.h"

template<>
InputParameters validParams<NSBndMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<NSBase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  MooseEnum bc_types("wall, far_field, symmetric, pressure_out, none", "none");  // 边界条件的类型，可以增加
  params.addRequiredParam<MooseEnum>("bc_type", bc_types, "边界条件");

  return params;
}

NSBndMaterial::NSBndMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		NSBase(name, parameters),
		_bc_type(getParam<MooseEnum>("bc_type")),
		_flux(declareProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable(declareProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable")),
		_flux_jacobi_grad_variable(declareProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable")),

		_penalty(declareProperty<std::vector<RealVectorValue> >("penalty")),
		_penalty_neighbor(declareProperty<std::vector<RealVectorValue> >("penalty_neighbor")),
		_penalty_jacobi_variable_ee(declareProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_ee")),
		_penalty_jacobi_variable_ne(declareProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_ne"))
{
	_n_equations = coupledComponents("variables");
	for (size_t eq = 0; eq < _n_equations; ++eq)
	{
		MooseVariable &val = *getVar("variables", eq);
		_ul.push_back(_is_implicit ? &val.sln() : &val.slnOld());
		_grad_ul.push_back(_is_implicit ? &val.gradSln(): &val.gradSlnOld());
	}
}

void NSBndMaterial::computeQpProperties()
{
	if(!_bnd)
	{
		mooseError("边界Material不在边界");
		return;
	}

	resizeQpProperty();

	Real ul[10], ur[10], ur_new[10];
	RealGradient dul[10], dur[10], dur_new[10];
	Real flux_new[10];
	RealVectorValue penalty_new[10], penalty_neighbor_new[10];

	computeQpLeftValue(ul, dul);
	computeQpRightValue(ur, dur, ul, dul);

	penaltyTerm(&_penalty[_qp][0], &_penalty_neighbor[_qp][0], ul, ur);
	fluxTerm(&_flux[_qp][0], ul, ur, dul, dur);

	for (int q = 0; q < _n_equations; ++q)
	{
		ul[q] += _ds;
		computeQpRightValue(ur_new, dur_new, ul, dul);
		penaltyTerm(penalty_new, penalty_neighbor_new, ul, ur_new);
		fluxTerm(flux_new, ul, ur_new, dul, dur_new);
		for (int p = 0; p < _n_equations; ++p)
		{
			_flux_jacobi_variable[_qp][p][q] = (flux_new[p] - _flux[_qp][p])/_ds;
			_penalty_jacobi_variable_ee[_qp][p][q] = (penalty_new[p] - _penalty[_qp][p])/_ds;
			_penalty_jacobi_variable_ne[_qp][p][q] = (penalty_neighbor_new[p] - _penalty_neighbor[_qp][p])/_ds;
		}
		ul[q] -= _ds;
	}

	for (int beta = 0; beta < 3; ++beta)
	for (int q = 0; q < _n_equations; ++q)
	{
		dul[q](beta) += _ds;
		computeQpRightValue(ur_new, dur_new, ul, dul);
		fluxTerm(flux_new, ul, ur_new, dul, dur_new);
		for (int p = 0; p < _n_equations; ++p)
		{
			_flux_jacobi_grad_variable[_qp][p][q](beta) = (flux_new[p] - _flux[_qp][p])/_ds;
		}
		dul[q](beta) -= _ds;
	}

}

void NSBndMaterial::computeQpLeftValue(Real* ul, RealGradient *dul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		ul[eq] = (*_ul[eq])[_qp];
		dul[eq] = (*_grad_ul[eq])[_qp];
	}
}

void NSBndMaterial::computeQpRightValue(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul)
{
	if(_bc_type == "wall")
	{
		wall(ur, dur, ul, dul);
		return;
	}
	if(_bc_type == "far_field")
	{
		farField(ur, dur, ul, dul);
		return;
	}
	if(_bc_type == "symmetric")
	{
		symmetric(ur, dur, ul, dul);
		return;
	}
	if(_bc_type == "pressure_out")
	{
		symmetric(ur, dur, ul, dul);
		return;
	}
	if(_bc_type == "none")
	{
		for (int eq = 0; eq < _n_equations; ++eq)
			ur[eq] = (*_ul[eq])[_qp];

		return;
	}

	mooseError(_bc_type<<"未定义的边界条件类型");
}

void NSBndMaterial::fluxTerm(Real *flux, Real* ul, Real* ur, RealGradient *dul, RealGradient *dur)
{
	RealVectorValue ifl[5], ifr[5], vfl[5], vfr[5];

	inviscousTerm(ifl, ul);
	inviscousTerm(ifr, ur);
	viscousTerm(vfl, ul, dul);
	viscousTerm(vfr, ur, dur);

	Real lam = (maxEigenValue(ul, _normals[_qp]) + maxEigenValue(ur, _normals[_qp]))/2.;
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux[eq] = 0.5*(ifl[eq] + ifr[eq] - (vfl[eq]+vfr[eq]))*_normals[_qp] + lam*(ul[eq] - ur[eq]);
	}
}

void NSBndMaterial::penaltyTerm(RealVectorValue* penalty, RealVectorValue* penalty_neighbor, Real* ul, Real* ur)
{
	RealGradient duh[10];
	for (int eq = 0; eq < _n_equations; ++eq)
		duh[eq] = (ul[eq]-ur[eq])/2.*_normals[_qp];

	viscousTerm(penalty, ul, duh);
	viscousTerm(penalty_neighbor, ur, duh);
}

void NSBndMaterial::wall(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		dur[eq] = dul[eq];

	const Point &normal = _normals[_qp];
	RealVectorValue momentum(ul[1], ul[2], ul[3]);
    Real vn = momentum*normal;
    Real pre = pressure(ul);
    Real twall = 1.;

    ur[0] = _gamma*_mach*_mach*pre/twall;
//    ur[0] = ul[0];
    ur[1] = 0.;
    ur[2] = 0.;
    ur[3] = 0.;
    ur[4] = pre/(_gamma-1) + 0.5*momentum.size_sq()/ur[0];
}

void NSBndMaterial::farField(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		dur[eq] = dul[eq];

	const Point &normal = _normals[_qp];

	Real rhoR, uR, vR, wR, pR;
	Real rhoL, uL, vL, wL, pL;
	Real cR, cL, cb;
	Real vnR, vnL, vnb;
	Real vel, s;
	Real Rp, Rm;

	uR = 1.0 * cos(_attack) * cos(_slide);
	vR = 1.0 * sin(_attack) * cos(_slide);
	wR = 1.0 * sin(_slide);
	rhoR = 1.0;

	pR = 1 / _gamma /_mach / _mach;
	cR = sqrt(fabs(_gamma * pR / rhoR));
	vnR = normal(0) * uR + normal(1) * vR + normal(2) * wR;

	rhoL = ul[0];
	uL = ul[1] / rhoL;
	vL = ul[2] / rhoL;
	wL = ul[3] / rhoL;
	vel = sqrt(uL * uL + vL * vL + wL * wL);
	pL = pressure(ul);
	cL = sqrt(fabs(_gamma * pL / rhoL));
	vnL =  normal(0) * uL + normal(1) * vL + normal(2) * wL;

	if (vel >= cL) {	//超声速
		if (vnL >= 0.0) //exit
		{
			ur[0] = ul[0];
			ur[1] = ul[1];
			ur[2] = ul[2];
			ur[3] = ul[3];
			ur[4] = ul[4];
		}
		else //inlet
		{
			ur[0] = rhoR;
			ur[1] = rhoR * uR;
			ur[2] = rhoR * vR;
			ur[3] = rhoR * wR;
			ur[4] = pR / (_gamma - 1) + 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);
		}
	}
	else
	{ 	//  亚声速
		if (vnL >= 0.0)
		{			//exit
			s = pL / pow(rhoL, _gamma);
			Rp = vnL + 2 * cL / (_gamma - 1);
			Rm = vnR - 2 * cR / (_gamma - 1);
			vnb = (Rp + Rm) / 2.0;
			cb = (Rp - Rm) * (_gamma - 1) / 4.0;

			ur[0] = pow((cb * cb) / (s * _gamma), 1.0 / (_gamma - 1));
			ur[1] = ur[0] * (uL + normal(0) * (vnb - vnL));
			ur[2] = ur[0] * (vL + normal(1) * (vnb - vnL));
			ur[3] = ur[0] * (wL + normal(2) * (vnb - vnL));
			ur[4] = cb * cb * ur[0] / _gamma / (_gamma - 1) + 0.5 * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]) / ur[0];
		}
		else
		{
			s = pR / pow(rhoR, _gamma);
			Rp = -vnR + 2.0 * cR / (_gamma - 1);
			Rm = -vnL - 2.0 * cL / (_gamma - 1);
			vnb = -(Rp + Rm) / 2.0;
			cb = (Rp - Rm) * (_gamma - 1) / 4.0;

			ur[0] = pow((cb * cb) / (s * _gamma), 1.0 / (_gamma - 1));
			ur[1] = ur[0] * (uR + normal(0) * (vnb - vnR));
			ur[2] = ur[0] * (vR + normal(1) * (vnb - vnR));
			ur[3] = ur[0] * (wR + normal(2) * (vnb - vnR));
			ur[4] = cb * cb * ur[0] / _gamma / (_gamma - 1) + 0.5 * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]) / ur[0];
		}
	}
}


void NSBndMaterial::symmetric(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		dur[eq] = dul[eq];

	const Point &normal = _normals[_qp];
	RealVectorValue momentum(ul[1], ul[2], ul[3]);
    Real vn = momentum*normal;
    Real pre = pressure(ul);

    ur[0] = ul[0];
    ur[1] = ul[1] - 2.0 * vn * normal(0);
    ur[2] = ul[2] - 2.0 * vn * normal(1);
    ur[3] = ul[3] - 2.0 * vn * normal(2);
    ur[4] = pre/(_gamma-1) + 0.5*momentum.size_sq()/ur[0];
}

void NSBndMaterial::pressureOut(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		dur[eq] = dul[eq];

	const Point &normal = _normals[_qp];
	RealVectorValue momentum(ul[1], ul[2], ul[3]);
    Real vn = momentum*normal;
    Real pre = 1 / _gamma /_mach / _mach;

    ur[0] = ul[0];
    ur[1] = ul[1];
    ur[2] = ul[2];
    ur[3] = ul[3];
    ur[4] = pre/(_gamma-1) + 0.5*momentum.size_sq()/ur[0];
}

void NSBndMaterial::resizeQpProperty()
{
	_flux[_qp].resize(_n_equations);
	_flux_jacobi_variable[_qp].resize(_n_equations);
	_flux_jacobi_grad_variable[_qp].resize(_n_equations);
	_penalty[_qp].resize(_n_equations);
	_penalty_neighbor[_qp].resize(_n_equations);
	_penalty_jacobi_variable_ee[_qp].resize(_n_equations);
	_penalty_jacobi_variable_ne[_qp].resize(_n_equations);
	for (int p = 0; p < _n_equations; ++p)
	{
		_flux_jacobi_variable[_qp][p].resize(_n_equations);
		_flux_jacobi_grad_variable[_qp][p].resize(_n_equations);
		_penalty_jacobi_variable_ee[_qp][p].resize(_n_equations);
		_penalty_jacobi_variable_ne[_qp][p].resize(_n_equations);
	}
}
