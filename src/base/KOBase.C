//
//#include "KOBase.h"
//
//template<>
//InputParameters validParams<KOBase>()
//{
//	InputParameters params = validParams<NSBase>();
//	return params;
//}
//
//KOBase::KOBase(const std::string& name, InputParameters parameters):
//		NSBase(name, parameters),
//		_sigma_k(0.5), _sigma_o(0.5), _beta_k(9./100), _beta_o(3./40), _alpha_o(5./9), _sigma_d(1./8), _c_lim(7./8),
//		_prandtl_turb(0.9),
//		_tu_infty(1e-09), _r_mu(1e-03)
//{
////	_omega_infty = _reynolds*1*_tu_infty/_r_mu;
//	_omega_infty = 1;
////	Real w_infty = 1000.;
//	_tu_infty = _omega_infty*_r_mu/_reynolds;
////	_r_mu = 1e-03;
//}
//
//void KOBase::inviscousTerm(RealVectorValue* inviscous_term, Real* uh)
//{
//	Real rho, p, h, k;
//	Real u, v, w;
//	rho = uh[0];
//	u = uh[1]/rho;
//	v = uh[2]/rho;
//	w = uh[3]/rho;
//	p = pressure(uh);
//	h = enthalpy(uh);
//	k = uh[5]/uh[0];
//
//	int component = 0;
//
//	component = 0;
//	inviscous_term[component](0) = uh[1];	// rhou
//	inviscous_term[component](1) = uh[2];	// rhov
//	inviscous_term[component](2) = uh[3];	// rhow
//
//	component = 1;
//	inviscous_term[component](0) = uh[1] * u + p + 2./3*uh[5];
//	inviscous_term[component](1) = uh[1] * v;
//	inviscous_term[component](2) = uh[1] * w;
//
//	component = 2;
//	inviscous_term[component](0) = uh[2] * u;
//	inviscous_term[component](1) = uh[2] * v + p + 2./3*uh[5];
//	inviscous_term[component](2) = uh[2] * w;
//
//	component = 3;
//	inviscous_term[component](0) = uh[3] * u;
//	inviscous_term[component](1) = uh[3] * v;
//	inviscous_term[component](2) = uh[3] * w + p + 2./3*uh[5];
//
//	component = 4;
//	inviscous_term[component](0) = rho * (h+2/3.*k) * u;
//	inviscous_term[component](1) = rho * (h+2/3.*k) * v;
//	inviscous_term[component](2) = rho * (h+2/3.*k) * w;
//
//	component = 5;
//	inviscous_term[component](0) = uh[5] * u;
//	inviscous_term[component](1) = uh[5] * v;
//	inviscous_term[component](2) = uh[5] * w;
////	inviscous_term[component].zero();
//
//	component = 6;
//	inviscous_term[component](0) = uh[6] * u;
//	inviscous_term[component](1) = uh[6] * v;
//	inviscous_term[component](2) = uh[6] * w;
////	inviscous_term[component].zero();
//}
//
//void KOBase::viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient* duh)
//{
//	Real rho = uh[0];
//	RealVectorValue velocity(uh[1]/rho, uh[2]/rho, uh[3]/rho);
//	RealGradient grad_rho(duh[0]);
//	RealTensor momentum_tensor(duh[1], duh[2], duh[3]);
//	RealTensor temp;
//	for (int i = 0; i < 3; ++i) {
//		for (int j = 0; j < 3; ++j)
//		{
//			temp(i,j) = velocity(i)*grad_rho(j);
//		}
//	}
//	RealTensor velocity_tensor = (momentum_tensor - temp)/rho;
//	RealTensor tau = velocity_tensor + velocity_tensor.transpose();
//	Real div = velocity_tensor(0,0) + velocity_tensor(1,1) + velocity_tensor(2,2);
//	Real lamdiv = 2./3. * div;
//	tau(0, 0) -= lamdiv; tau(1, 1) -= lamdiv; tau(2, 2) -= lamdiv;
//
//	Real mu = physicalViscosity(uh);
////	Real mu_turb = eddyViscosity(uh);
//	Real w_lim = log(7./8*sqrt(tau.size_sq()/_beta_k/2.));
//	Real k = uh[5]/uh[0];
//	Real w = uh[6]/uh[0];
//	w = std::max<Real>(w, w_lim);
//	Real mu_turb = _reynolds*uh[0]*k/exp(w);
//
//	tau *= (mu+mu_turb)/_reynolds;
//
//	RealVectorValue grad_enthalpy = (duh[4]-uh[4]/uh[0] * duh[0])/rho - velocity_tensor.transpose() * velocity;
//	grad_enthalpy *= (mu/_prandtl+mu_turb/_prandtl_turb)/_reynolds*(_gamma);
//	RealVectorValue grad_k = (duh[5] - uh[5]/uh[0]*duh[0])/uh[0]/_reynolds;
//	RealVectorValue grad_o = (duh[6] - uh[6]/uh[0]*duh[0])/uh[0]/_reynolds;
//
//	int component = 0;
//	viscous_term[component](0) = 0.;
//	viscous_term[component](1) = 0.;
//	viscous_term[component](2) = 0.;
//
//	component = 1;
//	viscous_term[component](0) = tau(0, 0);
//	viscous_term[component](1) = tau(0, 1);
//	viscous_term[component](2) = tau(0, 2);
//
//	component = 2;
//	viscous_term[component](0) = tau(1, 0);
//	viscous_term[component](1) = tau(1, 1);
//	viscous_term[component](2) = tau(1, 2);
//
//	component = 3;
//	viscous_term[component](0) = tau(2, 0);
//	viscous_term[component](1) = tau(2, 1);
//	viscous_term[component](2) = tau(2, 2);
//
//	component = 4;
//	RealVectorValue vel_tau = tau * velocity + grad_enthalpy ;
//	viscous_term[component](0) = vel_tau(0)+(mu + _sigma_k*mu_turb)*grad_k(0);
//	viscous_term[component](1) = vel_tau(1)+(mu + _sigma_k*mu_turb)*grad_k(1);
//	viscous_term[component](2) = vel_tau(2)+(mu + _sigma_k*mu_turb)*grad_k(2);
//
//	component = 5;
//
//	viscous_term[component](0) = (mu + _sigma_k*mu_turb)*grad_k(0);
//	viscous_term[component](1) = (mu + _sigma_k*mu_turb)*grad_k(1);
//	viscous_term[component](2) = (mu + _sigma_k*mu_turb)*grad_k(2);
////	viscous_term[component].zero();
//
//	component = 6;
//	viscous_term[component](0) = (mu + _sigma_o*mu_turb)*grad_o(0);
//	viscous_term[component](1) = (mu + _sigma_o*mu_turb)*grad_o(1);
//	viscous_term[component](2) = (mu + _sigma_o*mu_turb)*grad_o(2);
////	viscous_term[component].zero();
//}
//
//void KOBase::sourceTerm(Real* source_term, Real* uh, RealGradient* duh)
//{
//	Real rho = uh[0];
//	RealVectorValue velocity(uh[1]/rho, uh[2]/rho, uh[3]/rho);
//	RealGradient grad_rho(duh[0]);
//	RealTensor momentum_tensor(duh[1], duh[2], duh[3]);
//	RealTensor temp;
//	for (int i = 0; i < 3; ++i) {
//		for (int j = 0; j < 3; ++j)
//		{
//			temp(i,j) = velocity(i)*grad_rho(j);
//		}
//	}
//	RealTensor velocity_tensor = (momentum_tensor - temp)/rho;
//	RealTensor stress = (velocity_tensor + velocity_tensor.transpose())/2.;
////	Real mu_turb = eddyViscosity(uh);
//	Real mu = physicalViscosity(uh);
//	Real production = 2*stress.size_sq();
//
//	Real lamdiv = (velocity_tensor(0,0) + velocity_tensor(1,1) + velocity_tensor(2,2))/3.;
//	RealTensor stress1 = stress;
//	stress1(0, 0) -= lamdiv; stress1(1, 1) -= lamdiv; stress1(2, 2) -= lamdiv;
//	Real w_lim = log(7./8*sqrt(2*stress1.size_sq()/_beta_k));
//
//	Real k = uh[5]/uh[0];
//	Real w = uh[6]/uh[0];
//	k = std::max<Real>(k, 0);
//	w = std::max<Real>(w, w_lim);
//	Real mu_turb = _reynolds*uh[0]*k/exp(w);
//
//	RealVectorValue grad_o = (duh[6] - uh[6]/uh[0]*duh[0])/uh[0];
//	RealVectorValue grad_k = (duh[5] - uh[5]/uh[0]*duh[0])/uh[0];
//	Real grad_o_square = grad_o.size_sq();
//
//	int component = 0;
//	source_term[component] = 0;
//
//	component = 1;
//	source_term[component] = 0;
//
//	component = 2;
//	source_term[component] = 0;
//
//	component = 3;
//	source_term[component] = 0;
//
//	component = 4;
//	source_term[component] = 0.;
//
//	component = 5;
//	source_term[component] = mu_turb/_reynolds*production - _beta_k*rho*k*exp(w);
////	source_term[component] = 0.;
//
//	component = 6;
//	source_term[component] = _alpha_o*rho/exp(w)*production
//						   - _beta_o*rho*exp(w)
//						   + grad_o_square*(mu + _sigma_o*mu_turb)/_reynolds;
////	source_term[component] = 0.;
//
////	if(grad_k*grad_o > 0)
////		source_term[component] += _sigma_d*grad_k*grad_o/exp(w);
//}
//
//Real KOBase::eddyViscosity(Real* uh)
//{
//	Real k = uh[5]/uh[0];
//	Real w = uh[6]/uh[0];
//
//	k = std::max<Real>(k, 0);
//	k = std::min<Real>(k, 10000*_tu_infty);
//	return _reynolds*uh[0]*k/exp(w);
//}
//
//Real KOBase::pressure(Real* uh)
//{
////	return NSBase::pressure(uh);
//	return (_gamma-1)*(uh[4] - uh[5] - 0.5*(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0]);
//}
//
//int KOBase::equationIndex(const std::string &var_name)
//{
//	int eq = -1;
//	if(var_name == "rho")
//		eq = 0;
//	if(var_name == "momentum_x")
//		eq = 1;
//	if(var_name == "momentum_y")
//		eq = 2;
//	if(var_name == "momentum_z")
//		eq = 3;
//	if(var_name == "rhoe")
//		eq = 4;
//	if(var_name == "rhok")
//		eq = 5;
//	if(var_name == "rhoo")
//		eq = 6;
//
//	if(eq < 0)
//		mooseError("不可知的变量名");
//
//	return eq;
//}
