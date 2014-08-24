
#include "CFDBase.h"

template<>
InputParameters validParams<CFDBase>()
{
  InputParameters params = validParams<MooseObject>();
  params.addRequiredParam<Real>("mach",     "马赫数");
  params.addRequiredParam<Real>("reynolds", "雷诺数");

  params.addParam<Real>("gamma", 1.4, "比热比");
  params.addParam<Real>("prandtl", 0.72, "prandtl数");
  params.addParam<Real>("attack", 0, "攻角");
  params.addParam<Real>("slide", 0, "侧滑角");

  params.addParam<Real>("epsilon", -1, "对称项罚值，可以取-1, 0 , 1，分别对应SIP, IIP, NIP");
  params.addParam<Real>("sigma", 1, "通量罚值，默认值为6");
  params.addParam<Real>("ds", 1.490116119384766e-08, "微扰量");
  return params;
}

CFDBase::CFDBase(const std::string & name, InputParameters parameters)
{
	_gamma = (parameters.get<Real>("gamma"));
	_prandtl = (parameters.get<Real>("prandtl"));
	_mach = (parameters.get<Real>("mach"));
	_reynolds = (parameters.get<Real>("reynolds"));

	_attack = (parameters.get<Real>("attack"));
	_slide = (parameters.get<Real>("slide"));

	_epsilon = (parameters.get<Real>("epsilon"));
	_sigma = (parameters.get<Real>("sigma"));
	_ds = (parameters.get<Real>("ds"));
}

Real CFDBase::pressure(Real *uh)
{
	return (_gamma-1) * (uh[4] - 0.5*(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0]);  //
}
Real CFDBase::physicalViscosity(Real* uh)
{
	return 1.0;
}

Real CFDBase::enthalpy(Real *uh)
{
	Real p = pressure(uh);
	return (uh[4] + p)/uh[0];
}

Real CFDBase::temperature(Real* uh)
{
	Real p = pressure(uh);
	return _gamma*_mach*_mach*p/uh[0];
}

Real CFDBase::mach_local(Real* uh)
{
	Real vel = std::sqrt(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0];
	Real c = std::sqrt(temperature(uh))/_mach;
	return vel/c;
}

Real CFDBase::maxEigenValue(Real *uh, const Point &normal)
{
	RealVectorValue vel(uh[1]/uh[0], uh[2]/uh[0], uh[3]/uh[0]);
	Real c = std::sqrt(_gamma*pressure(uh)/uh[0]);
	return std::fabs(vel*normal)+c;
}

void CFDBase::inviscousTerm(RealVectorValue* inviscous_term, Real* uh)
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

//void CFDBase::fluxRiemann(Real* flux, Real* ul, Real* ur, Point& normal)
//{
//	RealVectorValue fl[5], fr[5];
//	inviscousTerm(fl, ul);
//	inviscousTerm(fr, ur);
//
//	Real rho, u, v, w, pre;
//	rho = (ul[0] + ur[0])/2.;
//	u = (ul[1] + ur [1])/rho;
//	v = (ul[2] + ur [2])/rho;
//	w = (ul[3] + ur [3])/rho;
//	pre = (pressure(ul) + pressure(ur))/2.;
//	Real lam = fabs(u*normal(0) + v * normal(1) + w * normal(2)) + sqrt(_gamma*pre/rho);
//	for (int eq = 0; eq < 5; ++eq)
//	{
//		flux[eq] = 0.5*((fl[eq] + fr[eq])*normal + lam*(ul[eq] - ur[eq]));
//	}
//
//}

void CFDBase::viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh)
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

void CFDBase::inviscousTerm(std::vector<RealVectorValue>& inviscous_term, Real* uh)
{
	inviscousTerm(&inviscous_term[0], uh);
}

//void CFDBase::liftOperator(Real *lift, Real *ul, Real *ur, Point &normal)
//{
//	Real uh[5];
//	RealGradient duh[5];
//	for(int eq = 0; eq < 5; ++eq)
//	{
//		uh[eq] = (ul[eq] + ur[eq])/2.;
//		duh[eq] = (ul[eq] - ur[eq])*normal;
//	}
//
//	RealVectorValue viscous_term[5];
//	viscousTerm(viscous_term, uh, duh);
//	for(int eq = 0; eq < 5; ++eq)
//	{
//		lift[eq] = viscous_term[eq]*normal;
//	}
////
////	Matrix5x5 G[3][3];
////	Real u, v, w, e, vel2;
////	u = uh[1]/uh[0];
////	v = uh[2]/uh[0];
////	w = uh[3]/uh[0];
////	e = uh[4]/uh[0];
////	vel2  = u*u+v*v+w*w;
////
////	G[0][0] <<
////			0, 0, 0, 0, 0,
////			-4./3*u, 4./3, 0, 0, 0,
////			-v,0,1,0,0,
////			-w,0,0,1,0,
////			-(vel2+1./3*u*u+_gamma/_prandtl*(e-vel2)), (4./3-_gamma/_prandtl)*u, (1-_gamma/_prandtl)*v, (1-_gamma/_prandtl)*w, _gamma/_prandtl;
////
////	G[0][1] <<
////			0,0,0,0,0,
////			2./3*v,0,-2./3,0,0,
////			-u,1,0,0,0,
////			0,0,0,0,0,
////			-1./3*u*v,v,-2./3*u,0,0;
////
////	G[0][2] <<
////			0,0,0,0,0,
////			2./3*w,0,0,-2./3,0,
////			0,0,0,0,0,
////			-u,-1,0,0,0,
////			-1./3*u*w,w,0,-2./3*u,0;
////
////	G[1][0] <<
////			0,0,0,0,0,
////			-v,0,1,0,0,
////			2./3*u,-2./3,0,0,0,
////			0,0,0,0,0,
////			-1./3*u*v,-2./3*v,u,0,0;
////
////	G[1][1] <<
////			0, 0, 0, 0, 0,
////			-u,1,0,0,0,
////			-4./3*v,0,4./3,0,0,
////			-w,0,0,1,0,
////			-(vel2+1./3*v*v+_gamma/_prandtl*(e-vel2)), (1-_gamma/_prandtl)*u, (4./3-_gamma/_prandtl)*v, (1-_gamma/_prandtl)*w, _gamma/_prandtl;
////
////	G[1][2] <<
////			0,0,0,0,0,
////			0,0,0,0,0,
////			2./3*w,0,0,-2./3,0,
////			-v,0,1,0,0,
////			-1./3*v*w,0,w,-2./3*v,0;
////
////	G[2][0] <<
////			0,0,0,0,0,
////			-w,0,0,1,0,
////			0,0,0,0,0,
////			2./3*u,-2./3,0,0,0,
////			-1./3*u*w,-2./3*w,0,u,0;
////
////	G[2][1] <<
////			0,0,0,0,0,
////			0,0,0,0,0,
////			-w,0,0,1,0,
////			2./3*v,0,-2./3,0,0,
////			-1./3*v*w,0,-2./3*w,v,0;
////
////	G[2][2] <<
////			0, 0, 0, 0, 0,
////			-u,1,0,0,0,
////			-v,0,1,0,0,
////			-4./3*w,0,0,4./3,0,
////			-(vel2+1./3*w*w+_gamma/_prandtl*(e-vel2)), (1-_gamma/_prandtl)*u, (1-_gamma/_prandtl)*v, (4./3-_gamma/_prandtl)*w, _gamma/_prandtl;
////
////	for (int i = 0; i < 5; ++i)
////	{
////		Real penalty = 0.;
////	for (int alpha = 0; alpha < 3; ++alpha)
////	{
////		for (int beta = 0; beta < 3; ++beta)
////		{
////			for (int j = 0; j < 5; ++j)
////			{
////				penalty += delta_uh[j]*normal(beta)*normal(alpha)*G[alpha][beta](i,j);
////			}
////		}
////	}
////	Real mu = physicalViscousity(uh);
////	lift[i] = penalty*mu/uh[0]/_reynolds;
////	}
//}
//
//Real CFDBase::penaltyOperator(Real* uh, Real* delta_uh, RealGradient& grad_phi, Point& normal, int i)
//{
//	Matrix5x5 G[3][3];
//	Real u, v, w, e, vel2;
//	u = uh[1]/uh[0];
//	v = uh[2]/uh[0];
//	w = uh[3]/uh[0];
//	e = uh[4]/uh[0];
//	vel2  = u*u+v*v+w*w;
//
//	G[0][0] <<
//			0, 0, 0, 0, 0,
//			-4./3*u, 4./3, 0, 0, 0,
//			-v,0,1,0,0,
//			-w,0,0,1,0,
//			-(vel2+1./3*u*u+_gamma/_prandtl*(e-vel2)), (4./3-_gamma/_prandtl)*u, (1-_gamma/_prandtl)*v, (1-_gamma/_prandtl)*w, _gamma/_prandtl;
//
//	G[0][1] <<
//			0,0,0,0,0,
//			2./3*v,0,-2./3,0,0,
//			-u,1,0,0,0,
//			0,0,0,0,0,
//			-1./3*u*v,v,-2./3*u,0,0;
//
//	G[0][2] <<
//			0,0,0,0,0,
//			2./3*w,0,0,-2./3,0,
//			0,0,0,0,0,
//			-u,-1,0,0,0,
//			-1./3*u*w,w,0,-2./3*u,0;
//
//	G[1][0] <<
//			0,0,0,0,0,
//			-v,0,1,0,0,
//			2./3*u,-2./3,0,0,0,
//			0,0,0,0,0,
//			-1./3*u*v,-2./3*v,u,0,0;
//
//	G[1][1] <<
//			0, 0, 0, 0, 0,
//			-u,1,0,0,0,
//			-4./3*v,0,4./3,0,0,
//			-w,0,0,1,0,
//			-(vel2+1./3*v*v+_gamma/_prandtl*(e-vel2)), (1-_gamma/_prandtl)*u, (4./3-_gamma/_prandtl)*v, (1-_gamma/_prandtl)*w, _gamma/_prandtl;
//
//	G[1][2] <<
//			0,0,0,0,0,
//			0,0,0,0,0,
//			2./3*w,0,0,-2./3,0,
//			-v,0,1,0,0,
//			-1./3*v*w,0,w,-2./3*v,0;
//
//	G[2][0] <<
//			0,0,0,0,0,
//			-w,0,0,1,0,
//			0,0,0,0,0,
//			2./3*u,-2./3,0,0,0,
//			-1./3*u*w,-2./3*w,0,u,0;
//
//	G[2][1] <<
//			0,0,0,0,0,
//			0,0,0,0,0,
//			-w,0,0,1,0,
//			2./3*v,0,-2./3,0,0,
//			-1./3*v*w,0,-2./3*w,v,0;
//
//	G[2][2] <<
//			0, 0, 0, 0, 0,
//			-u,1,0,0,0,
//			-v,0,1,0,0,
//			-4./3*w,0,0,4./3,0,
//			-(vel2+1./3*w*w+_gamma/_prandtl*(e-vel2)), (1-_gamma/_prandtl)*u, (1-_gamma/_prandtl)*v, (4./3-_gamma/_prandtl)*w, _gamma/_prandtl;
//
//	Real penalty = 0;
//	for (int alpha = 0; alpha < 3; ++alpha)
//	{
//		for (int beta = 0; beta < 3; ++beta)
//		{
//			for (int j = 0; j < 5; ++j)
//			{
//				penalty += delta_uh[j]*normal(beta)*grad_phi(alpha)*G[alpha][beta](i,j);
//			}
//		}
//	}
//	Real mu = physicalViscousity(uh);
//	return penalty*mu/uh[0]/_reynolds;
//}
