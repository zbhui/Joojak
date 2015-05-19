
#include "SAProblem.h"

using Eigen::Vector3d;

template<>
InputParameters validParams<SAProblem>()
{
  InputParameters params = validParams<NavierStokesProblem>();

  return params;
}

SAProblem::SAProblem(const std::string & name, InputParameters params) :
	NavierStokesProblem(name, params),
	_cb1(0.1355), _cb2(0.622), _sigma_sa(2./3), _kappa(0.41),
	_cw2(0.3), _cw3(2.0), _cv1(7.1), _cv2(0.7), _cv3(0.9),
	_ct1(1.0), _ct2(2.0), _ct3(1.2), _ct4(0.5),
	_prandtl_turb(0.9),
	_nu_infty(3.0)
{
	_cw1 = _cb1/_kappa/_kappa + (1+_cb2)/_sigma_sa;
	_cw3_pow6 = _cw3*_cw3*_cw3*_cw3*_cw3*_cw3;

	_n_equations = 6;
}

void SAProblem::inviscousTerm(RealVectorValue* inviscous_term, Real* uh)
{
	NavierStokesProblem::inviscousTerm(inviscous_term, uh);
	Real rho, u, v, w;
	rho = uh[0];
	u = uh[1]/rho;
	v = uh[2]/rho;
	w = uh[3]/rho;

	int component = 5;
	inviscous_term[component](0) = uh[5] * u;
	inviscous_term[component](1) = uh[5] * v;
	inviscous_term[component](2) = uh[5] * w;
}

void SAProblem::viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh)
{
	Real rho(uh[0]);
	RealVectorValue velocity(uh[1]/rho, uh[2]/rho, uh[3]/rho);
	RealGradient grad_rho(duh[0]);
	RealTensor momentum_tensor(duh[1], duh[2], duh[3]);
	RealTensor temp;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j)
		{
			temp(i,j) = velocity(i)*grad_rho(j);
		}
	}
	RealTensor velocity_tensor((momentum_tensor - temp)/rho);
	RealTensor tau(velocity_tensor + velocity_tensor.transpose());
	Real lamdiv = 2./3. * velocity_tensor.tr();
	tau(0, 0) -= lamdiv; tau(1, 1) -= lamdiv; tau(2, 2) -= lamdiv;

	Real mu = physicalViscosity(uh);
	Real X = uh[5]/mu;
	Real psi(X);
	if (X < 10)
		psi = 0.05*log(1+exp(20*X));
	else
		psi = X;

	Real fv1 = psi*psi*psi/(psi*psi*psi+_cv1*_cv1*_cv1);
	Real mu_turb = std::max<Real>(0., uh[5]*fv1);

	tau *= (mu+mu_turb)/_reynolds;

	RealVectorValue grad_enthalpy = (duh[4]-uh[4]/uh[0] * duh[0])/uh[0] - velocity_tensor.transpose() * velocity;
	RealVectorValue grad_nu = (duh[5] - uh[5]/uh[0]*duh[0])/uh[0];

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
	RealVectorValue vel_tau(tau * velocity + (mu/_prandtl+mu_turb/_prandtl_turb)/_reynolds*(_gamma)*grad_enthalpy);
	viscous_term[component](0) = vel_tau(0);
	viscous_term[component](1) = vel_tau(1);
	viscous_term[component](2) = vel_tau(2);

	component = 5;
	viscous_term[component](0) = mu*(1+psi)/_sigma_sa*grad_nu(0)/_reynolds;
	viscous_term[component](1) = mu*(1+psi)/_sigma_sa*grad_nu(1)/_reynolds;
	viscous_term[component](2) = mu*(1+psi)/_sigma_sa*grad_nu(2)/_reynolds;

}

void SAProblem::viscousTermAdiabatic(RealVectorValue* viscous_term, Real* uh, RealGradient *duh)
{
	Real rho(uh[0]);
	RealVectorValue velocity(uh[1]/rho, uh[2]/rho, uh[3]/rho);
	RealGradient grad_rho(duh[0]);
	RealTensor momentum_tensor(duh[1], duh[2], duh[3]);
	RealTensor temp;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j)
		{
			temp(i,j) = velocity(i)*grad_rho(j);
		}
	}
	RealTensor velocity_tensor((momentum_tensor - temp)/rho);
	RealTensor tau(velocity_tensor + velocity_tensor.transpose());
	Real lamdiv = 2./3. * velocity_tensor.tr();
	tau(0, 0) -= lamdiv; tau(1, 1) -= lamdiv; tau(2, 2) -= lamdiv;

	Real mu = physicalViscosity(uh);
	Real X = uh[5]/mu;
	Real psi(X);
	if (X < 10)
		psi = 0.05*log(1+exp(20*X));
	else
		psi = X;

	Real fv1 = psi*psi*psi/(psi*psi*psi+_cv1*_cv1*_cv1);
	Real mu_turb = std::max<Real>(0., uh[5]*fv1);

	tau *= (mu+mu_turb)/_reynolds;

	RealVectorValue grad_enthalpy = (duh[4]-uh[4]/uh[0] * duh[0])/uh[0] - velocity_tensor.transpose() * velocity;
	RealVectorValue grad_nu = (duh[5] - uh[5]/uh[0]*duh[0])/uh[0];

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
	RealVectorValue vel_tau(tau * velocity);
	viscous_term[component](0) = vel_tau(0);
	viscous_term[component](1) = vel_tau(1);
	viscous_term[component](2) = vel_tau(2);

	component = 5;
	viscous_term[component](0) = mu*(1+psi)/_sigma_sa*grad_nu(0)/_reynolds;
	viscous_term[component](1) = mu*(1+psi)/_sigma_sa*grad_nu(1)/_reynolds;
	viscous_term[component](2) = mu*(1+psi)/_sigma_sa*grad_nu(2)/_reynolds;

}

void SAProblem::sourceTerm(Real* source_term, Real* uh, RealGradient* duh)
{
	Real rho(uh[0]);
	Real d(uh[6]);     //uh[6] 为壁面距离
	RealVectorValue velocity(uh[1]/rho, uh[2]/rho, uh[3]/rho);
	RealGradient grad_rho(duh[0]);
	RealTensor momentum_tensor(duh[1], duh[2], duh[3]);
	RealTensor temp;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j)
		{
			temp(i,j) = velocity(i)*grad_rho(j);
		}
	}
	RealTensor velocity_tensor((momentum_tensor - temp)/rho);
	RealTensor tau(velocity_tensor + velocity_tensor.transpose());
	Real lamdiv = 2./3. * velocity_tensor.tr();
	tau(0, 0) -= lamdiv; tau(1, 1) -= lamdiv; tau(2, 2) -= lamdiv;

	Real mu = physicalViscosity(uh);
	Real X = uh[5]/mu;
	Real psi(X);
	if (X < 10)
		psi = 0.05*log(1+exp(20*X));
	else
		psi = X;

	Real fv1 = psi*psi*psi/(psi*psi*psi+_cv1*_cv1*_cv1);
	Real mu_turb = std::max<Real>(0., uh[5]*fv1);

	tau *= (mu+mu_turb)/_reynolds;

	RealVectorValue grad_enthalpy = (duh[4]-uh[4]/uh[0] * duh[0])/uh[0] - velocity_tensor.transpose() * velocity;
	RealVectorValue grad_nu = (duh[5] - uh[5]/uh[0]*duh[0])/uh[0];


	RealTensor omega((velocity_tensor - velocity_tensor.transpose())/2.);
	Real vorticity = sqrt(2*omega.size_sq())+1E-04;

	Real fv2, s_title, s_hat;
	fv2 = 1-psi/(1+psi*fv1);
	s_hat = mu*psi/(_reynolds*rho*_kappa*_kappa*d*d)*fv2;
	if(s_hat >= -_cv2*vorticity)
		s_title = vorticity+s_hat;
	else
		s_title = vorticity+(_cv2*_cv2*vorticity+_cv3*s_hat)*vorticity/((_cv3-2*_cv2)*vorticity-s_hat);

	Real r = std::min<Real>((s_hat)/((s_title*fv2)), 10);
	//	Real r = s_hat/s_title/fv2;
	//	std::cout << r  <<std::endl;
	Real r6 = r*r*r*r*r*r;
	Real g = r+_cw2*(r6-r);
	Real g6 = g*g*g*g*g*g;
	Real fw = g*pow((1+_cw3_pow6)/(g6+_cw3_pow6),1./6);

	int component = 0;
	source_term[component] = 0;

	component = 1;
	source_term[component] = 0;

	component = 2;
	source_term[component] = 0;

	component = 3;
	source_term[component] = 0;

	component = 4;
	source_term[component] = 0.;

	component = 5;
	source_term[component] = _cb1*s_title*mu*psi
						   + _cb2/_sigma_sa/_reynolds*rho*grad_nu.size_sq()
						   - _cw1/_reynolds/rho*fw*(mu*mu*psi*psi/d/d);
}

void SAProblem::fluxRiemann(Real* flux, Real* ul, Real* ur, Point &normal)
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

void SAProblem::boundaryCondition(Real *ur, Real *ul, Point &normal, std::string bc_type)
{
	if(bc_type == "isothermal_wall")
	{
		isothermalWall(ur, ul, normal);
		return;
	}
	if(bc_type == "adiabatic_wall")
	{
		adiabaticWall(ur, ul, normal);
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
//	if(_bc_type == "symmetric")
//	{
//		symmetric(ur, dur, ul, dul);
//		return;
//	}
//	if(_bc_type == "pressure_out")
//	{
//		symmetric(ur, dur, ul, dul);
//		return;
//	}
//	if(_bc_type == "none")
//	{
//		for (int eq = 0; eq < _n_equations; ++eq)
//			ur[eq] = (*_ul[eq])[_qp];
//
//		return;
//	}

	mooseError( bc_type << "未定义的边界条件类型");
}

void SAProblem::computeBoundaryFlux(Real* flux, RealVectorValue* lift, Real* ul, RealGradient* dul, Point& normal, Real penalty, std::string bc_type)
{
	Real ur[10];
	RealGradient dur[10];
	RealVectorValue ifl[10], ifr[10], vfl[10], vfr[10];

	for (int eq = 0; eq < _n_equations; ++eq)
		dur[eq] = dul[eq];

	boundaryCondition(ur, ul, normal, bc_type);
	computeLift(lift, ul, ur, normal);

	if(bc_type == "adiabatic_wall")
	{
		viscousTermAdiabatic(vfl, ul, dul);
		viscousTermAdiabatic(vfr, ur, dur);
	}
	else
	{
		viscousTerm(vfl, ul, dul);
		viscousTerm(vfr, ur, dur);
	}

	fluxRiemann(flux, ul, ur, normal);

	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux[eq] -= 0.5*((vfl[eq]+vfr[eq])-penalty*lift[eq])*normal;
	}
}

void SAProblem::isothermalWall(Real *ur,  Real *ul, Point &normal)
{
    NavierStokesProblem::isothermalWall(ur, ul, normal);
    ur[5] = 0;
}

void SAProblem::adiabaticWall(Real *ur,  Real *ul, Point &normal)
{
    NavierStokesProblem::adiabaticWall(ur, ul, normal);
    ur[5] = 0;
}

void SAProblem::symmetric(Real *ur,  Real *ul, Point &normal)
{
    NavierStokesProblem::symmetric(ur, ul, normal);
    ur[5] = ul[5];
}
void SAProblem::farField(Real *ur, Real *ul, Point &normal)
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
	Real rho = ul[0];
	uL = ul[1] / rhoL;
	vL = ul[2] / rhoL;
	wL = ul[3] / rhoL;
	vel = sqrt(uL * uL + vL * vL + wL * wL);
	pL = pressure(ul);
	cL = sqrt(fabs(_gamma * pL / rhoL));
	vnL =  normal(0) * uL + normal(1) * vL + normal(2) * wL;

	if(lam[1] < 0)  // 入口
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
		ur[5] = rho*_nu_infty;
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
		ur[5] = ul[5];
	}
}

Real SAProblem::eddyViscosity(Real* uh)
{
	if(uh[5] < 0)
		return 0.;

	Real mu = physicalViscosity(uh);
	Real X = uh[5]/mu;
	Real psi;
	if (X <= 10)
		psi = 0.05*log(1+exp(20*X));
	else
		psi = X;

	Real fv1 = psi*psi*psi/(psi*psi*psi+_cv1*_cv1*_cv1);
	return uh[5]*fv1;
}

Real SAProblem::computeAuxValue(std::string var_name, Real* uh)
{
	if(var_name == "pressure")
		return pressure(uh);
	if(var_name == "mach")
		return localMach(uh);
	if(var_name == "velocity_x")
		return uh[1]/uh[0];
	if(var_name == "velocity_y")
		return uh[2]/uh[0];
	if(var_name == "velocity_z")
		return uh[3]/uh[0];
	if(var_name == "eddy_viscosity")
		return eddyViscosity(uh);

	mooseError("未知的辅助变量");
	return 0;
}
