
#include "NavierStokesProblem.h"

#include "MooseApp.h"
using namespace Eigen;

template<>
InputParameters validParams<NavierStokesProblem>()
{
  InputParameters params = validParams<CLawProblem>();
  params.addParam<Real>("mach",  0.1, "马赫数");
  params.addParam<Real>("gamma", 1.4, "比热比");
  params.addParam<Real>("reynolds", 1, "雷诺数");
  params.addParam<Real>("prandtl", 0.72, "prandtl数");
  params.addParam<Real>("attack", 0., "攻角");
  params.addParam<Real>("sideslip", 0., "侧滑角");
  params.addParam<Real>("pitch", 0., "俯仰角");
  params.addParam<Real>("yaw", 180., "偏航角");
  params.addParam<Real>("roll", -90., "滚转角");

  params.addParam<Real>("ref_length", 1, "参考长度");
  params.addParam<Real>("ref_area", 1, "参考面积");

  return params;
}

NavierStokesProblem::NavierStokesProblem(const std::string & name, InputParameters params) :
	CLawProblem(name, params),
	_mach(getParam<Real>("mach")),
	_gamma(getParam<Real>("gamma")),
	_reynolds(getParam<Real>("reynolds")),
	_prandtl(getParam<Real>("prandtl")),

	_attack(getParam<Real>("attack")*libMesh::pi/180),
	_sideslip(getParam<Real>("sideslip")*libMesh::pi/180),
	_pitch(getParam<Real>("pitch")*libMesh::pi/180),
	_yaw(getParam<Real>("yaw")*libMesh::pi/180),
	_roll(getParam<Real>("roll")*libMesh::pi/180),

	_ref_length(getParam<Real>("ref_length")),
	_ref_area(getParam<Real>("ref_area"))
{
	_n_equations = 5;
}

Real NavierStokesProblem::pressure(Real *uh)
{
	return (_gamma-1) * (uh[4] - 0.5*(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0]);  //
}

Real NavierStokesProblem::pressureInfity()
{
	return 1.0/_gamma/_mach/_mach;
}
Real NavierStokesProblem::enthalpy(Real *uh)
{
	Real p = pressure(uh);
	return (uh[4] + p)/uh[0];
}

Real NavierStokesProblem::temperature(Real* uh)
{
	Real p = pressure(uh);
	return _gamma*_mach*_mach*p/uh[0];
}

Real NavierStokesProblem::mach_local(Real* uh)
{
	Real vel = std::sqrt(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0];
	Real c = std::sqrt(temperature(uh))/_mach;
	return vel/c;
}

Real NavierStokesProblem::maxEigenValue(Real *uh, const Point &normal)
{
	RealVectorValue vel(uh[1]/uh[0], uh[2]/uh[0], uh[3]/uh[0]);
	Real c = std::sqrt(_gamma*pressure(uh)/uh[0]);
	return std::fabs(vel*normal)+c;
}

void NavierStokesProblem::eigenValue(Real *lam, Real *uh, const Point &normal)
{
	RealVectorValue vel(uh[1]/uh[0], uh[2]/uh[0], uh[3]/uh[0]);
	Real c = std::sqrt(_gamma*pressure(uh)/uh[0]);
	lam[0] = vel*normal-c;
	lam[1] = vel*normal;
	lam[2] = vel*normal+c;
}

void NavierStokesProblem::inviscousTerm(RealVectorValue* inviscous_term, Real* uh)
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

void NavierStokesProblem::viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh)
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
	Real mu = physicalViscosity(uh);
	tau *= mu/_reynolds;

	RealVectorValue grad_enthalpy = (duh[4]-uh[4]/uh[0] * duh[0])/rho - velocity_tensor.transpose() * velocity;
	grad_enthalpy *= (mu/_reynolds)*(_gamma/_prandtl);

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

void NavierStokesProblem::inviscousTerm(std::vector<RealVectorValue>& inviscous_term, Real* uh)
{
	inviscousTerm(&inviscous_term[0], uh);
}

Eigen::Quaterniond NavierStokesProblem::bodyFromWind()
{
	Quaterniond q_attack(AngleAxisd(_attack, Vector3d::UnitY()));
	Quaterniond q_sideslip (AngleAxisd(pi-_sideslip,  Vector3d::UnitZ()));
	return q_attack*q_sideslip;
}

Eigen::Quaterniond NavierStokesProblem::earthFromBody()
{
	Quaterniond q_pitch(AngleAxisd(_pitch, Vector3d::UnitY()));
	Quaterniond q_yaw(AngleAxisd(_yaw, Vector3d::UnitZ()));
	Quaterniond q_roll(AngleAxisd(_roll, Vector3d::UnitX()));
	return q_roll*q_pitch*q_yaw;
}

Eigen::Quaterniond NavierStokesProblem::earthFromWind()
{
	return earthFromBody()*bodyFromWind();
}

Real NavierStokesProblem::physicalViscosity(Real* uh)
{
	return 1.0;
}

int NavierStokesProblem::equationIndex(const std::string &var_name)
{
	int eq = -1;
	if(var_name == "rho")
		eq = 0;
	if(var_name == "momentum_x")
		eq = 1;
	if(var_name == "momentum_y")
		eq = 2;
	if(var_name == "momentum_z")
		eq = 3;
	if(var_name == "rhoe")
		eq = 4;

	if(eq < 0)
		mooseError("不可知的变量名");

	return eq;
}

void NavierStokesProblem::fluxRiemann(Real* flux, Real* ul, Real* ur, const Point &normal)
{
	RealVectorValue ifl[5], ifr[5], vfl[5], vfr[5];

	inviscousTerm(ifl, ul);
	inviscousTerm(ifr, ur);

	Real lam = 1;//(maxEigenValue(ul, _normals[_qp]) + maxEigenValue(ur, _normals[_qp]))/2.;
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux[eq] = 0.5*(ifl[eq] + ifr[eq])*normal + lam*(ul[eq] - ur[eq]);
	}
}

void NavierStokesProblem::boundaryCondition(Real *ur, Real *ul, Point &normal, MooseEnum _bc_type)
{
//	switch (_bc_type)
//	{
//		case "isothermal_wall":
//			isothermalWall(ur, ul, normal);
//			break;
//		case "adiabatic_wall":
//			isothermalWall(ur, ul, normal);
//			break;
//		case "far_field":
//			isothermalWall(ur, ul, normal);
//			break;
//		case "symmetric":
//			isothermalWall(ur, ul, normal);
//			break;
//		case "pressure_out":
//			isothermalWall(ur, ul, normal);
//			break;
//		case "none":
//			isothermalWall(ur, ul, normal);
//			break;
//		default:
//			mooseError(_bc_type<<"未定义的边界条件类型");
//			break;
//	}
}

void NavierStokesProblem::isothermalWall(Real *ur,  Real *ul, Point &normal)
{
    Real twall = 1.;
    Real pre = pressure(ul);

    ur[0] = ul[0];
    ur[1] = 0.;
    ur[2] = 0.;
    ur[3] = 0.;
    ur[4] = ul[0]*twall/_gamma/(_gamma-1)/_mach/_mach;
}

//void NavierStokesProblem::adiabaticWall(Real* ur, RealGradient* dur, Real* ul, RealGradient* dul)
//{
//	for (int eq = 0; eq < _n_equations; ++eq)
//		dur[eq] = dul[eq];
//
//    ur[0] = ul[0];
//    ur[1] = 0.;
//    ur[2] = 0.;
//    ur[3] = 0.;
//    ur[4] = ul[4];
//}
//
//void NavierStokesProblem::farFieldRiemann(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul)
//{
//	for (int eq = 0; eq < _n_equations; ++eq)
//		dur[eq] = dul[eq];
//
//	const Point &normal = _normals[_qp];
//
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
//
//void NavierStokesProblem::farField(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul)
//{
//	for (int eq = 0; eq < _n_equations; ++eq)
//		dur[eq] = dul[eq];
//
//	const Point &normal = _normals[_qp];
//
//	Real rhoR, uR, vR, wR, pR;
//	Real rhoL, uL, vL, wL, pL;
//	Real cR, cL, cb;
//	Real vnR, vnL, vnb;
//	Real vel, s;
//	Real Rp, Rm;
//
//	Vector3d vel_inf = earthFromWind()*Vector3d::UnitX();
//	if(_current_elem->dim() == 2)
//		vel_inf(2) = 0.;
//	uR = vel_inf(0);
//	vR = vel_inf(1);
//	wR = vel_inf(2);
//
//	Real lam[3];
//	eigenValue(lam, ul, normal);
//
//	rhoR = 1.0;
//
//	pR = 1 / _gamma /_mach / _mach;
//	cR = sqrt(fabs(_gamma * pR / rhoR));
//	vnR = normal(0) * uR + normal(1) * vR + normal(2) * wR;
//
//	rhoL = ul[0];
//	uL = ul[1] / rhoL;
//	vL = ul[2] / rhoL;
//	wL = ul[3] / rhoL;
//	vel = sqrt(uL * uL + vL * vL + wL * wL);
//	pL = pressure(ul);
//	cL = sqrt(fabs(_gamma * pL / rhoL));
//	vnL =  normal(0) * uL + normal(1) * vL + normal(2) * wL;
//
//	if(lam[1] < 0)  // 入口
//	{
//		if (vel > cL) //超音速
//		{
//			ur[0] = rhoR;
//			ur[1] = rhoR * uR;
//			ur[2] = rhoR * vR;
//			ur[3] = rhoR * wR;
//			ur[4] = pR / (_gamma - 1) + 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);
//		}
//		else	//亚音速
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
//	else  //出口
//	{
//		if (vel > cL) //超音速
//		{
//			ur[0] = ul[0];
//			ur[1] = ul[1];
//			ur[2] = ul[2];
//			ur[3] = ul[3];
//			ur[4] = ul[4];
//		}
//		else	//亚音速
//		{
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
//	}
//}
//
//
//void NavierStokesProblem::symmetric(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul)
//{
//	for (int eq = 0; eq < _n_equations; ++eq)
//		dur[eq] = dul[eq];
//
//	const Point &normal = _normals[_qp];
//	RealVectorValue momentum(ul[1], ul[2], ul[3]);
//    Real vn = momentum*normal;
//    Real pre = pressure(ul);
//
//    ur[0] = ul[0];
//    ur[1] = ul[1] - 2.0 * vn * normal(0);
//    ur[2] = ul[2] - 2.0 * vn * normal(1);
//    ur[3] = ul[3] - 2.0 * vn * normal(2);
//    ur[4] = ul[4];
////    ur[4] = pre/(_gamma-1) + 0.5*momentum.size_sq()/ur[0];
//}
//
//void NavierStokesProblem::pressureOut(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul)
//{
//	for (int eq = 0; eq < _n_equations; ++eq)
//		dur[eq] = dul[eq];
//
//	Vector3d vel_left(ul[1]/ul[0], ul[2]/ul[0], ul[3]/ul[0]);
//	Real p_inf = 1/_gamma/_mach/_mach;
//	ur[0] = ul[0];
//	ur[1] = ul[1];
//	ur[2] = ul[2];
//	ur[3] = ul[3];
//	ur[4] = p_inf/(_gamma - 1) + 0.5*ur[0]*vel_left.squaredNorm();
//}
