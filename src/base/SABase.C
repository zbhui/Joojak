
#include "SABase.h"

template<>
InputParameters validParams<SABase>()
{
	InputParameters params = validParams<SABase>();
	return params;
}

SABase::SABase(const std::string& name, InputParameters parameters):
		NSBase(name, parameters),
		_cb1(0.1355), _cb2(0.622), _sigma_sa(2./3), _kappa(0.41),
		_cw2(0.3), _cw3(2.0), _cv1(7.1), _cv2(0.7), _cv3(0.9),
		_ct1(1.0), _ct2(2.0), _ct3(1.2), _ct4(0.5),
		_prandtl_turb(0.9),
		_nu_infty(0.01)
{
	_cw1 = _cb1/_kappa/_kappa + (1+_cb2)/_sigma_sa;
}

void SABase::inviscousTerm(RealVectorValue* inviscous_term, Real* uh)
{
	Real rho, p, h, k;
	Real u, v, w;
	rho = uh[0];
	u = uh[1]/rho;
	v = uh[2]/rho;
	w = uh[3]/rho;
	p = pressure(uh);
	h = enthalpy(uh);
	k = uh[5]/uh[0];

	int component = 0;

	component = 0;
	inviscous_term[component](0) = uh[1];
	inviscous_term[component](1) = uh[2];
	inviscous_term[component](2) = uh[3];

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
	inviscous_term[component](0) = rho*h * u;
	inviscous_term[component](1) = rho*h * v;
	inviscous_term[component](2) = rho*h * w;

	component = 5;
	inviscous_term[component](0) = uh[5] * u;
	inviscous_term[component](1) = uh[5] * v;
	inviscous_term[component](2) = uh[5] * w;
}

void SABase::viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient* duh)
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
	Real X = uh[5]/mu;
	Real psi;
	if (X <= 10)
		psi = 0.05*log(1+exp(20*X));
	else
		psi = X;

	Real fv1 = pow(psi,3)/(pow(psi,3)+pow(_cv1,3));
	Real mu_turb = mu*psi*fv1;

	tau *= (mu+mu_turb)/_reynolds;

	RealVectorValue grad_enthalpy = (duh[4]-uh[4]/uh[0] * duh[0])/rho - velocity_tensor.transpose() * velocity;
	grad_enthalpy *= (mu/_prandtl+mu_turb/_prandtl_turb)/_reynolds*(_gamma);
	RealVectorValue grad_nu = (duh[5] - uh[5]/uh[0]*duh[0])/uh[0]/_reynolds;

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
	viscous_term[component](0) = mu*(1+psi)*grad_nu(0);
	viscous_term[component](1) = mu*(1+psi)*grad_nu(1);
	viscous_term[component](2) = mu*(1+psi)*grad_nu(2);
}

void SABase::sourceTerm(Real* source_term, Real* uh, RealGradient* duh, Real d)
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
	RealTensor stress = (velocity_tensor + velocity_tensor.transpose())/2.;
	RealTensor omega =  (velocity_tensor - velocity_tensor.transpose())/2.;
	Real vorticity = sqrt(2*omega.size_sq());

	Real mu = physicalViscosity(uh);
	Real X = uh[5]/mu;
	Real psi;
	if (X <= 10)
		psi = 0.05*log(1+exp(20*X));
	else
		psi = X;

	Real fv1, fv2, s_title, s_hat;
	fv1 = pow(psi,3)/(pow(psi,3)+pow(_cv1,3));
	fv2 = 1-psi/(1+psi*fv1);
	s_hat = mu*psi/(_reynolds*rho*_kappa*_kappa*d*d)*fv2;
	if(s_hat >= -_cv2*vorticity)
		s_title = vorticity+s_hat;
	else
		s_title = vorticity+(_cv2*_cv2*vorticity+_cv3*s_hat)*vorticity/((_cv3-2*_cv2)*vorticity-s_hat);

	Real r = std::min<Real>(s_hat/s_title/fv2, 10);
	Real g = r+_cw2*(pow(r,6)-r);
	Real fw = g*(pow((1+pow(_cw3, 6))/((pow(g,6)+pow(_cw3, 6))),1./6));

	RealVectorValue grad_nu = (duh[5] - uh[5]/uh[0]*duh[0])/uh[0];

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
						   + 1./_sigma_sa/_reynolds*_cb2*rho*grad_nu.size_sq()
						   - _cw1/_reynolds/rho*fw*(mu*mu*psi*psi/d/d);
}

Real SABase::eddyViscosity(Real* uh)
{
	Real mu = physicalViscosity(uh);
	Real X = uh[5]/mu;
	Real psi;
	if (X <= 10)
		psi = 0.05*log(1+exp(20*X));
	else
		psi = X;

	Real fv1 = pow(psi,3)/(pow(psi,3)+pow(_cv1,3));
	return mu*psi*fv1;
}

int SABase::equationIndex(const std::string &var_name)
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
	if(var_name == "rhon")
		eq = 5;

	if(eq < 0)
		mooseError("不可知的变量名");

	return eq;
}
