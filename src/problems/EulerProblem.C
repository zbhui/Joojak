
#include "EulerProblem.h"
#include "CLawCellMaterial.h"
#include "CLawFaceMaterial.h"
#include "CLawBoundaryMaterial.h"

using Eigen::Vector3d;

template<>
InputParameters validParams<EulerProblem>()
{
  InputParameters params = validParams<CFDProblem>();

  return params;
}

EulerProblem::EulerProblem(const std::string & name, InputParameters params) :
	CFDProblem(name, params)
{
}

void EulerProblem::computeCellFlux(RealGradient* flux, Real* source, Real* uh, RealGradient* duh)
{
	RealVectorValue inv_term[10], vis_term[10], source_term[10], artificial_term[10];
	inviscousTerm(inv_term, uh);
	viscousTerm(vis_term, uh, duh);
	sourceTerm(source, uh, duh);
	artificialViscous(artificial_term, uh, duh);
	for (int eq = 0; eq < _n_equations; ++eq)
		flux[eq] = inv_term[eq] - vis_term[eq] - artificial_term[eq];
}

void EulerProblem::computeFaceFlux(Real* flux, RealVectorValue* lift, Real* ul, Real* ur, RealGradient* dul, RealGradient* dur, Point& normal, Real penalty)
{
	RealVectorValue ifl[10], ifr[10], vfl[10], vfr[10];
	computeLift(lift, ul, ur, normal);

	viscousTerm(vfl, ul, dul);
	viscousTerm(vfr, ur, dur);

	fluxRiemann(flux, ul, ur, normal);

	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux[eq] -= 0.5*((vfl[eq]+vfr[eq])-penalty*lift[eq])*normal;
	}

}
void EulerProblem::artificialViscous(RealVectorValue* artificial_viscous, Real* uh, RealGradient* duh)
{
	RealVectorValue vel(uh[1], uh[2], uh[3]);
	vel /= uh[0];
	Real c = std::sqrt(temperature(uh))/_mach;
    Real indicator = uh[5] ;
	Real viscosity = indicator;//* duh[0].size()/uh[0]*0.01;
//	Real viscosity = 0.00;//4*(vel.size() + c)*indicator;

	int component = 0;
	artificial_viscous[component] = viscosity*duh[component];

	component = 1;
	artificial_viscous[component] = viscosity*duh[component];

	component = 2;
	artificial_viscous[component] = viscosity*duh[component];

	component = 3;
	artificial_viscous[component] = viscosity*duh[component];

	component = 4;
	artificial_viscous[component] = viscosity*duh[component];
}

void EulerProblem::inviscousTerm(RealVectorValue* inviscous_term, Real* uh)
{
	Real rho, p, h;
	Real u, v, w;
	rho = uh[0];
	u = uh[1]/rho;
	v = uh[2]/rho;
	w = uh[3]/rho;
	p = pressure(uh);
	h = enthalpy(uh);

	int component = 0;

	component = 0;
	inviscous_term[component](0) = uh[1];	// rhou
	inviscous_term[component](1) = uh[2];	// rhov
	inviscous_term[component](2) = uh[3];	// rhow

	component = 1;
	inviscous_term[component](0) = uh[1] * u + p;
	inviscous_term[component](1) = uh[1] * v;
	inviscous_term[component](2) = uh[1] * w;

	component = 2;
	inviscous_term[component](0) = uh[2] * u;
	inviscous_term[component](1) = uh[2] * v + p;
	inviscous_term[component](2) = uh[2] * w;

	component = 3;
	inviscous_term[component](0) = uh[3] * u;
	inviscous_term[component](1) = uh[3] * v;
	inviscous_term[component](2) = uh[3] * w + p;

	component = 4;
	inviscous_term[component](0) = rho * h * u;
	inviscous_term[component](1) = rho * h * v;
	inviscous_term[component](2) = rho * h * w;
}

void EulerProblem::viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient* duh)
{

	Real viscosity = 0.000;
	int component = 0;
	viscous_term[component](0) = 0.;
	viscous_term[component](1) = 0.;
	viscous_term[component](2) = 0.;
	viscous_term[component] = viscosity*duh[component];

	component = 1;
	viscous_term[component](0) = 0.;
	viscous_term[component](1) = 0.;
	viscous_term[component](2) = 0.;
	viscous_term[component] = viscosity*duh[component];

	component = 2;
	viscous_term[component](0) = 0.;
	viscous_term[component](1) = 0.;
	viscous_term[component](2) = 0.;
	viscous_term[component] = viscosity*duh[component];

	component = 3;
	viscous_term[component](0) = 0.;
	viscous_term[component](1) = 0.;
	viscous_term[component](2) = 0.;
	viscous_term[component] = viscosity*duh[component];

	component = 4;
	viscous_term[component](0) = 0.;
	viscous_term[component](1) = 0.;
	viscous_term[component](2) = 0.;
	viscous_term[component] = viscosity*duh[component];
}

void EulerProblem::fluxRiemann(Real* flux, Real* ul, Real* ur, Point &normal)
{
	RealVectorValue ifl[5], ifr[5], vfl[5], vfr[5];
	Real uh[5];

	inviscousTerm(ifl, ul);
	inviscousTerm(ifr, ur);

	for (int eq = 0; eq < _n_equations; ++eq)
		uh[eq] = (ul[eq]+ur[eq])/2;

	for (int eq = 0; eq < _n_equations; ++eq)
		flux[eq] = 0.5*(ifl[eq] + ifr[eq])*normal + maxEigenValue(uh, normal)*(ul[eq] - ur[eq]);
}

void EulerProblem::boundaryCondition(Real *ur, Real *ul, Point &normal, std::string bc_type)
{
	if(bc_type == "wall")
	{
		wall(ur, ul, normal);
		return;
	}
	if(bc_type == "far_field")
	{
		farField(ur, ul, normal);
		return;
	}
	if(bc_type == "symmetric")
	{
		symmetric(ur, ul, normal);
		return;
	}
	mooseError( bc_type << "未定义的边界条件类型");
}

void EulerProblem::wall(Real *ur,  Real *ul, Point &normal)
{
	RealVectorValue momentum(ul[1], ul[2], ul[3]);
    Real vn = momentum*normal;
    Real pre = pressure(ul);

    ur[0] = ul[0];
    ur[1] = ul[1] - 2.0 * vn * normal(0);
    ur[2] = ul[2] - 2.0 * vn * normal(1);
    ur[3] = ul[3] - 2.0 * vn * normal(2);
    ur[4] = pre/(_gamma-1) + 0.5*momentum.size_sq()/ur[0];
}

void EulerProblem::symmetric(Real *ur,  Real *ul, Point &normal)
{
	wall(ur, ul, normal);
}

void EulerProblem::farField(Real *ur,  Real *ul, Point &normal)
{
	Real rhoR, uR, vR, wR, pR;
	Real rhoL, uL, vL, wL, pL;
	Real cR, cL, cb;
	Real vnR, vnL, vnb;
	Real vel, s;
	Real Rp, Rm;

	Vector3d vel_inf = _attitude.earthFromWind()*Vector3d::UnitX();
	if(_mesh.dimension() == 2)
		vel_inf(2) = 0.;

	uR = vel_inf(0);
	vR = vel_inf(1);
	wR = vel_inf(2);

	Real lam[3];
	eigenValue(lam, ul, normal);

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

	if(vnR <= 0)  // 入口
	{
		if (vel > cL) //超音速
		{
			ur[0] = rhoR;
			ur[1] = rhoR * uR;
			ur[2] = rhoR * vR;
			ur[3] = rhoR * wR;
			ur[4] = pR / (_gamma - 1) + 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);
		}
		else	//亚音速
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
	else  //出口
	{
		if (vel > cL) //超音速
		{
			ur[0] = ul[0];
			ur[1] = ul[1];
			ur[2] = ul[2];
			ur[3] = ul[3];
			ur[4] = ul[4];
		}
		else	//亚音速
		{
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
	}
//	DependValue left, right;
//	Real rhoR, uR, vR, wR, pR;
//	Real rhoL, uL, vL, wL, pL;
//	Real cR, cL, cb;
//	Real vnR, vnL, vnb;
//	Real vel, s;
//	Real Rp, Rm;
//
//	Vector3d velocity = _attitude.earthFromWind()*Vector3d::UnitX();
//	uR = velocity(0);
//	vR = velocity(1);
//	wR = velocity(2);
//	rhoR = 1.0;
//
//	right.rho = 1;
//	right.p = 1 / _gamma /_mach / _mach;
//	left.vel = RealVectorValue(velocity(0), velocity(1), velocity(2));
//	right.c = _mach * 1;
//
//    left.rho = ul[0];
//    left.vel = RealVectorValue(ul[1], ul[2], ul[3])/ul[0];
//    left.p = pressure(ul);
//    left.c = sqrt(fabs(_gamma * left.p / left.rho));
//
//	Real vn_left = left.vel*normal;
//	Real vn_right = right.vel*normal;
//	vel = left.vel.size();
//
//	if (vel >= cL) {	//超声速
//		if (vnL >= 0.0) //exit
//		{
//			ur[0] = ul[0];
//			ur[1] = ul[1];
//			ur[2] = ul[2];
//			ur[3] = ul[3];
//			ur[4] = ul[4];
//		}
//		else //inlet
//		{
//			ur[0] = right.rho;
//			ur[1] = right.rho*right.vel(0);
//			ur[2] = right.rho*right.vel(1);
//			ur[3] = right.rho*right.vel(2);
//			ur[4] = right.p/(_gamma - 1) + 0.5*right.rho*right.vel.size();
//		}
//	}
//	else
//	{ 	//  亚声速
//		if (vnL >= 0.0)
//		{			//exit
//			s = pL / pow(rhoL, _gamma);
//			Rp = vnL + 2 * cL / (_gamma - 1);
//			Rm = vnR - 2 * cR / (_gamma - 1);
//			vnb = (Rp + Rm) / 2.0;
//			cb = (Rp - Rm) * (_gamma - 1) / 4.0;
//
//			ur[0] = pow((cb * cb) / (s * _gamma), 1.0 / (_gamma - 1));
//			ur[1] = ur[0] * (uL + normal(0) * (vnb - vnL));
//			ur[2] = ur[0] * (vL + normal(1) * (vnb - vnL));
//			ur[3] = ur[0] * (wL + normal(2) * (vnb - vnL));
//			ur[4] = cb * cb * ur[0] / _gamma / (_gamma - 1) + 0.5 * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]) / ur[0];
//		}
//		else
//		{
//			s = pR / pow(rhoR, _gamma);
//			Rp = -vnR + 2.0 * cR / (_gamma - 1);
//			Rm = -vnL - 2.0 * cL / (_gamma - 1);
//			vnb = -(Rp + Rm) / 2.0;
//			cb = (Rp - Rm) * (_gamma - 1) / 4.0;
//
//			ur[0] = pow((cb * cb) / (s * _gamma), 1.0 / (_gamma - 1));
//			ur[1] = ur[0] * (uR + normal(0) * (vnb - vnR));
//			ur[2] = ur[0] * (vR + normal(1) * (vnb - vnR));
//			ur[3] = ur[0] * (wR + normal(2) * (vnb - vnR));
//			ur[4] = cb * cb * ur[0] / _gamma / (_gamma - 1) + 0.5 * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]) / ur[0];
//		}
//	}
}

//void EulerProblem::farFieldRiemann(Real *ur,  Real *ul, Point &normal)
//{
//	Vector3d vel_inf = earthFromWind()*Vector3d::UnitX();
//	if(_current_elem->dim() == 2)
//		vel_inf(2) = 0.;
//
//	Real rho_inf = 1.;
//	Real p_inf = 1/_gamma/_mach/_mach;
//	Real pl = pressure(ul);
//	Vector3d vel_left(ul[1]/ul[0], ul[2]/ul[0], ul[3]/ul[0]);
////	Real vel = vel_left.norm();
//	Real vel = vel_inf.norm();
//	Real cl = sqrt(fabs(_gamma * pl /ul[0]));
////	Real vn = vel_left(0)*normal(0)+vel_left(1)*normal(1)+vel_left(2)*normal(2);
//	Real vn = vel_inf(0)*normal(0)+vel_inf(1)*normal(1)+vel_inf(2)*normal(2);
//
//	if(vn < 0)  // 入口
//	{
//		if (vel > cl) //超音速
//		{
//			ur[0] = rho_inf;
//			ur[1] = rho_inf*vel_inf(0);
//			ur[2] = rho_inf*vel_inf(1);
//			ur[3] = rho_inf*vel_inf(2);
//			ur[4] = p_inf/(_gamma - 1) + 0.5 * rho_inf*vel_inf.squaredNorm();
//		}
//		else	//亚音速
//		{
//			ur[0] = rho_inf;
//			ur[1] = rho_inf*vel_inf(0);
//			ur[2] = rho_inf*vel_inf(1);
//			ur[3] = rho_inf*vel_inf(2);
//			ur[4] = pl/(_gamma - 1) + 0.5 * rho_inf*vel_inf.squaredNorm();
//		}
//	}
//	else  //出口
//	{
//		if (vel > cl) //超音速
//		{
//			ur[0] = ul[0];
//			ur[1] = ul[1];
//			ur[2] = ul[2];
//			ur[3] = ul[3];
//			ur[4] = ul[4];
//		}
//		else	//亚音速
//		{
//			ur[0] = ul[0];
//			ur[1] = ul[1];
//			ur[2] = ul[2];
//			ur[3] = ul[3];
//			ur[4] = p_inf/(_gamma - 1) + 0.5*ur[0]*vel_left.squaredNorm();
//		}
//	}
//}
Real EulerProblem::physicalViscosity(Real* uh)
{
	return 0;
}


void EulerProblem::updateDependValue(DependValue& denpend_value, Real *uh)
{
	denpend_value.p = pressure(uh);
}
