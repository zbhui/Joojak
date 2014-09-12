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

#include "EulerBndMaterial.h"

template<>
InputParameters validParams<EulerBndMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<EulerBase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  MooseEnum bc_types("wall, far_field, symmetric, pressure_out, none", "none");  // 边界条件的类型，可以增加
  params.addRequiredParam<MooseEnum>("bc_type", bc_types, "边界条件");

  return params;
}

EulerBndMaterial::EulerBndMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		EulerBase(name, parameters),
		_bc_type(getParam<MooseEnum>("bc_type")),
		_flux(declareProperty<std::vector<Real> >("flux")),
		_jacobi_variable(declareProperty<std::vector<std::vector<Real> > >("bnd_jacobi_variable"))
{
	_n_equations = coupledComponents("variables");
	for (size_t eq = 0; eq < _n_equations; ++eq)
	{
		MooseVariable &val = *getVar("variables", eq);
		_ul.push_back(_is_implicit ? &val.sln() : &val.slnOld());
	}
}

void EulerBndMaterial::computeQpProperties()
{
	if(!_bnd)
	{
		mooseError("边界Material不在边界");
		return;
	}

	resizeQpProperty();

	Real ul[10], ur[10], ul_new[10], ur_new[10];
	Real flux_new[10];

	computeQpLeftValue(ul);
	computeQpRightValue(ur);
	fluxRiemann(&_flux[_qp][0], ul, ur);

	Matrix5x5 jacobi_variable_en, jacobi_variable_ur_ul;
	for (int q = 0; q < _n_equations; ++q)
	{
//		ul[q] += _ds;
		(*_ul[q])[_qp] += _ds;
		computeQpLeftValue(ul);
		computeQpRightValue(ur_new);
		fluxRiemann(flux_new, ul, ur_new);
		for (int p = 0; p < _n_equations; ++p)
			_jacobi_variable[_qp][p][q] = (flux_new[p] - _flux[_qp][p])/_ds;

//		ul[q] -= _ds;
		(*_ul[q])[_qp] -= _ds;

//		ur[q] += _ds;
//		fluxRiemann(flux_new, ul, ur);
//		for (int p = 0; p < _n_equations; ++p)
//		{
//			Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
//			jacobi_variable_en(p,q) = tmp;
//		}
//		ur[q] -= _ds;
//
//		(*_ul[q])[_qp] += _ds;
//		computeQpRightValue(ur_new);
//		for (int p = 0; p < _n_equations; ++p)
//		{
//			Real tmp = (ur_new[p] - ur[p])/_ds;
//			jacobi_variable_ur_ul(p,q) = tmp;
//		}
//		(*_ul[q])[_qp] -= _ds;
	}
//
////	std::cout << jacobi_variable_en*jacobi_variable_ur_ul <<std::endl;
//	for (int q = 0; q < _n_equations; ++q)
//	{
//		for (int p = 0; p < _n_equations; ++p)
//		{
//			_jacobi_variable[_qp][p][q] += (jacobi_variable_en*jacobi_variable_ur_ul)(p,q);
//		}
//	}
}

void EulerBndMaterial::computeQpLeftValue(Real* ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_ul[eq])[_qp];
}

void EulerBndMaterial::computeQpRightValue(Real* ur)
{
	if(_bc_type == "wall")
	{
		wall(ur);
		return;
	}
	if(_bc_type == "far_field")
	{
		farField(ur);
		return;
	}
	if(_bc_type == "symmetric")
	{
		symmetric(ur);
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

void EulerBndMaterial::fluxRiemann(Real *flux, Real* ul, Real* ur)
{
	const Point &normal = _normals[_qp];
	RealVectorValue fl[10], fr[10];
	inviscousTerm(fl, ul);
	inviscousTerm(fr, ur);

	Real lam = (maxEigenValue(ul, _normals[_qp]) + maxEigenValue(ur, _normals[_qp]))/2.;
//	lam = 1.;
	for (int eq = 0; eq < _n_equations; ++eq)
		flux[eq] = 0.5*(fl[eq] + fr[eq])*normal + lam*(ul[eq] - ur[eq]);
}

void EulerBndMaterial::wall(Real* ur)
{
	Real ul[5];
	computeQpLeftValue(ul);

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

void EulerBndMaterial::farField(Real* ur)
{
	Real ul[5];
	computeQpLeftValue(ul);

	const Point &normal = _normals[_qp];

	Real rhoR, uR, vR, wR, pR;
	Real rhoL, uL, vL, wL, pL;
	Real cR, cL, cb;
	Real vnR, vnL, vnb;
	Real vel, s;
	Real Rp, Rm;

	uR = 1.0 * cos(_attack) * cos(_sideslip);
	vR = 1.0 * sin(_attack) * cos(_sideslip);
	wR = 1.0 * sin(_sideslip);
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

void EulerBndMaterial::symmetric(Real* ur)
{
	wall(ur);
}

void EulerBndMaterial::resizeQpProperty()
{
	_flux[_qp].resize(_n_equations);

	_jacobi_variable[_qp].resize(_n_equations);
	for (int p = 0; p < _n_equations; ++p)
		_jacobi_variable[_qp][p].resize(_n_equations);
}
