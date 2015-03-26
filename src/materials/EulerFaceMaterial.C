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
  params += validParams<EulerBase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

EulerFaceMaterial::EulerFaceMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		EulerBase(name, parameters),
		_flux(declareProperty<std::vector<Real> >("flux")),
		_jacobi_variable_ee(declareProperty<std::vector<std::vector<Real> > >("face_jacobi_ee")),
		_jacobi_variable_en(declareProperty<std::vector<std::vector<Real> > >("face_jacobi_en")),
		_jacobi_variable_ne(declareProperty<std::vector<std::vector<Real> > >("face_jacobi_ne")),
		_jacobi_variable_nn(declareProperty<std::vector<std::vector<Real> > >("face_jacobi_nn"))
{
	_n_equations = coupledComponents("variables");

	if(_bnd && _neighbor)
	{
		for (int eq = 0; eq < _n_equations; ++eq)
		{
			MooseVariable &val = *getVar("variables", eq);
			_ul.push_back(_is_implicit ? &val.sln() : &val.slnOld());
			_ur.push_back(_is_implicit ? &val.slnNeighbor() : &val.slnOldNeighbor());
			_grad_ul.push_back(_is_implicit ? &val.gradSln(): &val.gradSlnOld());
			_grad_ur.push_back(_is_implicit ? &val.gradSlnNeighbor(): &val.gradSlnOldNeighbor());
		}
	}
}

void EulerFaceMaterial::computeQpProperties()
{
	if(_bnd && _neighbor)
	{
		_flux[_qp].resize(_n_equations);

		Real ul[10], ur[10];
		RealGradient dul[10], dur[10];
		Real flux_new[10];

		computeQpLeftValue(ul);
		computeQpRightValue(ur);
		computeQpLeftGradValue(dul);
		computeQpRightGradValue(dur);
		fluxTerm(&_flux[_qp][0], ul, ur, dul, dur);

		_jacobi_variable_ee[_qp].resize(_n_equations);
		_jacobi_variable_en[_qp].resize(_n_equations);
		_jacobi_variable_ne[_qp].resize(_n_equations);
		_jacobi_variable_nn[_qp].resize(_n_equations);
		for (int p = 0; p < _n_equations; ++p)
		{
			_jacobi_variable_ee[_qp][p].resize(_n_equations);
			_jacobi_variable_en[_qp][p].resize(_n_equations);
			_jacobi_variable_ne[_qp][p].resize(_n_equations);
			_jacobi_variable_nn[_qp][p].resize(_n_equations);
		}

		for (int q = 0; q < _n_equations; ++q)
		{
			ul[q] += _ds;
			fluxTerm(flux_new, ul, ur, dul, dur);
			for (int p = 0; p < _n_equations; ++p)
			{
				Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
				_jacobi_variable_ee[_qp][p][q] = tmp;
				_jacobi_variable_ne[_qp][p][q] = -tmp;
			}
			ul[q] -= _ds;

			ur[q] += _ds;
			fluxTerm(flux_new, ul, ur, dul, dur);
			for (int p = 0; p < _n_equations; ++p)
			{
				Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
				_jacobi_variable_en[_qp][p][q] = tmp;
				_jacobi_variable_nn[_qp][p][q] = -tmp;
			}
			ur[q] -= _ds;
		}
	}
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

void EulerFaceMaterial::computeQpLeftGradValue(RealGradient *ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_grad_ul[eq])[_qp];
}

void EulerFaceMaterial::computeQpRightGradValue(RealGradient *ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_grad_ur[eq])[_qp];
}

void EulerFaceMaterial::artificalViscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient* duh)
{
	Real rho = uh[0];
	RealVectorValue velocity(uh[1]/rho, uh[2]/rho, uh[3]/rho);
	RealGradient grad_rho(duh[0]);
	RealTensor momentum_tensor(duh[1], duh[2], duh[3]);
	RealTensor temp;
	for (int alpha = 0; alpha < 3; ++alpha) {
		for (int beta = 0; beta < 3; ++beta)
		{
			temp(alpha,beta) = velocity(alpha)*grad_rho(beta);
		}
	}
	RealTensor velocity_tensor = (momentum_tensor - temp)/rho;
	RealTensor tau = velocity_tensor + velocity_tensor.transpose();
	Real div = velocity_tensor(0,0) + velocity_tensor(1,1) + velocity_tensor(2,2);
	Real lamdiv = 2./3. * div;
	tau(0, 0) -= lamdiv; tau(1, 1) -= lamdiv; tau(2, 2) -= lamdiv;
	Real mu = artificalViscous(uh);
	tau *= mu;

	Real prandl = 0.72;
	RealVectorValue grad_enthalpy = (duh[4]-uh[4]/uh[0] * duh[0])/rho - velocity_tensor.transpose() * velocity;
	grad_enthalpy *= (mu)*(_gamma/prandl);

	int component = 0;
	viscous_term[component](0) = 0.;
	viscous_term[component](1) = 0.;
	viscous_term[component](2) = 0.;

	component = 1;
	viscous_term[component](0) = tau(0, 0);
	viscous_term[component](1) = tau(0, 1);
	viscous_term[component](2) = tau(0, 2);

	component = 2;
	viscous_term[component](0) = tau(1, 0);
	viscous_term[component](1) = tau(1, 1);
	viscous_term[component](2) = tau(1, 2);

	component = 3;
	viscous_term[component](0) = tau(2, 0);
	viscous_term[component](1) = tau(2, 1);
	viscous_term[component](2) = tau(2, 2);

	component = 4;
	RealVectorValue vel_tau = tau * velocity + grad_enthalpy ;
	viscous_term[component](0) = vel_tau(0);
	viscous_term[component](1) = vel_tau(1);
	viscous_term[component](2) = vel_tau(2);
}

Real EulerFaceMaterial::artificalViscous(Real* uuh)
{
	const double h = _assembly.elemVolume()/_assembly.sideElemVolume();
    Real c = 0.2;

	Real ul[10], ur[10], uh[10];
	computeQpLeftValue(ul);
	computeQpRightValue(ur);
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		uh[eq] = 0.5*(ul[eq]+ur[eq]);
	}

	RealVectorValue fl[10], fr[10];
	inviscousTerm(fl, ul);
	inviscousTerm(fr, ur);

	Real sensitive(0);
	Real dp_du[10];
	dp_du[0] = uh[4]/uh[0];
	dp_du[1] = -uh[1]/uh[0];
	dp_du[2] = -uh[2]/uh[0];
	dp_du[3] = -uh[3]/uh[0];
	dp_du[4] = 1.;

	Real pre = pressure(uh);
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		sensitive += fabs((fl[eq] - fr[eq])*_normals[_qp]*dp_du[eq]);
	}

	return c*h*h*_grad_ul[0]->size()*sensitive/pre;
}

void EulerFaceMaterial::fluxTerm(Real *flux, Real* ul, Real* ur, RealGradient *dul, RealGradient *dur)
{
	RealVectorValue ifl[5], ifr[5], vfl[5], vfr[5];

	inviscousTerm(ifl, ul);
	inviscousTerm(ifr, ur);
	artificalViscousTerm(vfl, ul, dul);
	artificalViscousTerm(vfr, ur, dur);

	Real lam = (maxEigenValue(ul, _normals[_qp]) + maxEigenValue(ur, _normals[_qp]))/2.;
	for (int eq = 0; eq < _n_equations; ++eq)
	{
//		flux[eq] = 0.5*(ifl[eq] + ifr[eq] - (vfl[eq]+vfr[eq]))*_normals[_qp] + lam*(ul[eq] - ur[eq]);
		flux[eq] = 0.5*(ifl[eq] + ifr[eq])*_normals[_qp] + lam*(ul[eq] - ur[eq]);
	}
}
