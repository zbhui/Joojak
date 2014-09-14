
#include "KOBase.h"

template<>
InputParameters validParams<KOBase>()
{
	InputParameters params = validParams<NSBase>();
	return params;
}

KOBase::KOBase(const std::string& name, InputParameters parameters):
		NSBase(name, parameters),
		_sigma_k(0.5), _sigma_o(0.5), _beta_k(0.09), _beta_o(0.0075), _alpha_o(5./9),
		_prandtl_turb(0.9),
		_tu_infty(0.005), _r_mu(1e-05)
{
}

void KOBase::inviscousTerm(RealVectorValue* inviscous_term, Real* uh)
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
	inviscous_term[component](0) = uh[1] * u + p + 2./3*uh[5];
	inviscous_term[component](1) = uh[1] * v;
	inviscous_term[component](2) = uh[1] * w;

	component = 2;
	inviscous_term[component](0) = uh[2] * u;
	inviscous_term[component](1) = uh[2] * v + p + 2./3*uh[5];
	inviscous_term[component](2) = uh[2] * w;

	component = 3;
	inviscous_term[component](0) = uh[3] * u;
	inviscous_term[component](1) = uh[3] * v;
	inviscous_term[component](2) = uh[3] * w + p + 2./3*uh[5];

	component = 4;
	inviscous_term[component](0) = rho * h * u;
	inviscous_term[component](1) = rho * h * v;
	inviscous_term[component](2) = rho * h * w;

	component = 5;
	inviscous_term[component](0) = uh[5] * u;
	inviscous_term[component](1) = uh[5] * v;
	inviscous_term[component](2) = uh[5] * w;

	component = 6;
	inviscous_term[component](0) = uh[6] * u;
	inviscous_term[component](1) = uh[6] * v;
	inviscous_term[component](2) = uh[6] * w;
}

void KOBase::viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient* duh)
{
	Real rho = uh[0];
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
	RealTensor velocity_tensor = (momentum_tensor - temp)/rho;
	RealTensor tau = velocity_tensor + velocity_tensor.transpose();
	Real div = velocity_tensor(0,0) + velocity_tensor(1,1) + velocity_tensor(2,2);
	Real lamdiv = 2./3. * div;
	tau(0, 0) -= lamdiv; tau(1, 1) -= lamdiv; tau(2, 2) -= lamdiv;
	Real mu = physicalViscosity(uh);
	Real mu_turb = eddyViscosity(uh);
	tau *= (mu+mu_turb)/_reynolds;

	RealVectorValue grad_enthalpy = (duh[4]-uh[4]/uh[0] * duh[0])/rho - velocity_tensor.transpose() * velocity;
	grad_enthalpy *= (mu/_prandtl+mu_turb/_prandtl_turb)/_reynolds*(_gamma);


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

	component = 5;
	RealVectorValue grad_k = (duh[5] - uh[5]/uh[0]*duh[0])/uh[0]/_reynolds;
	viscous_term[component](0) = (mu + _sigma_k*mu_turb)*grad_k(0);
	viscous_term[component](1) = (mu + _sigma_k*mu_turb)*grad_k(1);
	viscous_term[component](2) = (mu + _sigma_k*mu_turb)*grad_k(2);

	component = 6;
	RealVectorValue grad_o = (duh[6] - uh[6]/uh[0]*duh[0])/uh[0]/_reynolds;
	viscous_term[component](0) = (mu + _sigma_o*mu_turb)*grad_o(0);
	viscous_term[component](1) = (mu + _sigma_o*mu_turb)*grad_o(1);
	viscous_term[component](2) = (mu + _sigma_o*mu_turb)*grad_o(2);
}

void KOBase::sourceTerm(Real* source_term, Real* uh, RealGradient* duh)
{
	Real rho = uh[0];
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
	RealTensor velocity_tensor = (momentum_tensor - temp)/rho;
	RealTensor sigma = velocity_tensor + velocity_tensor.transpose();
	Real div = velocity_tensor(0,0) + velocity_tensor(1,1) + velocity_tensor(2,2);
	Real lamdiv = 2./3. * div;
	sigma(0, 0) -= lamdiv; sigma(1, 1) -= lamdiv; sigma(2, 2) -= lamdiv;
	Real mu_turb = eddyViscosity(uh);
	Real mu = physicalViscosity(uh);
	sigma *= mu_turb/_reynolds;
	sigma(0, 0) -= 2./3*uh[5]; sigma(1, 1) -= 2./3*uh[5]; sigma(2, 2) -= 2./3*uh[5];

	Real source = 0;
	for (int a = 0; a < LIBMESH_DIM; ++a)
	{
		for (int b = 0; b < LIBMESH_DIM; ++b)
		{
			source += sigma(a,b)*velocity_tensor(a,b);
		}
	}

	Real k = uh[5]/uh[0];
	Real w = uh[6]/uh[0];
	RealVectorValue grad_o = (duh[6] - uh[6]/uh[0]*duh[0])/uh[0];
	Real grad_o_square = grad_o.size_sq();

	int component = 0;
	source_term[component] = 0;

	component = 1;
	source_term[component] = 0;

	component = 2;
	source_term[component] = 0;

	component = 3;
	source_term[component] = 0;

	component = 4;
	source_term[component] = -source + _beta_k*rho*k*std::exp(w);

	component = 5;
	source_term[component] = source - _beta_k*rho*k*std::exp(w);

	component = 6;
	source_term[component] = (mu + _sigma_o*mu_turb)*grad_o_square/_reynolds + _alpha_o/k*source - _beta_o*rho*std::exp(w);
}

Real KOBase::eddyViscosity(Real* uh)
{
	return _reynolds * uh[5]/std::exp(uh[6]/uh[0]);
}


