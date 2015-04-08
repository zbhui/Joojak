
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
