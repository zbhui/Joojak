
#include "NSBase.h"

template<>
InputParameters validParams<NSBase>()
{
  InputParameters params = validParams<EulerBase>();
  params.addParam<Real>("reynolds", 10, "雷诺数");
  params.addParam<Real>("prandtl", 0.72, "prandtl数");
  params.addParam<Real>("epsilon", 1, "对称项罚值，可以取1, 0 , -1，分别对应SIP, IIP, NIP");
  params.addParam<Real>("sigma", 1, "通量罚值，默认值为6");
  return params;
}

NSBase::NSBase(const std::string & name, InputParameters parameters):
		EulerBase(name, parameters)
{
	_prandtl = (parameters.get<Real>("prandtl"));
	_reynolds = (parameters.get<Real>("reynolds"));

	_epsilon = (parameters.get<Real>("epsilon"));
	_sigma = (parameters.get<Real>("sigma"));
}


Real NSBase::physicalViscosity(Real* uh)
{
	return 1.0;
}

void NSBase::viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh)
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

void NSBase::stressTerm(RealTensorValue &tau, Real* uh, RealGradient* duh)
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
	tau = velocity_tensor + velocity_tensor.transpose();
	Real div = velocity_tensor(0,0) + velocity_tensor(1,1) + velocity_tensor(2,2);
	Real lamdiv = 2./3. * div;
	tau(0, 0) -= lamdiv; tau(1, 1) -= lamdiv; tau(2, 2) -= lamdiv;
	Real mu = physicalViscosity(uh);
	tau *= mu/_reynolds;
}

void NSBase::stressTerm(Matrix3d &tau, Real* uh, RealGradient* duh)
{
	RealTensorValue tmp;
	stressTerm(tmp, uh, duh);
	for (int alpha = 0; alpha < 3; ++alpha) {
		for (int beta = 0; beta < 3; ++beta)
		{
			tau(alpha,beta) = tmp(alpha, beta);
		}
	}
}

int NSBase::equationIndex(const std::string &var_name)
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
