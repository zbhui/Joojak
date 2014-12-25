
#include "EulerBase.h"

template<>
InputParameters validParams<EulerBase>()
{
  InputParameters params = validParams<MooseObject>();
  params.addParam<Real>("mach",  0.2, "马赫数");
  params.addParam<Real>("gamma", 1.4, "比热比");
  params.addParam<Real>("attack", 0., "攻角");
  params.addParam<Real>("sideslip", 0., "侧滑角");
  params.addParam<Real>("pitch", 0., "俯仰角");
  params.addParam<Real>("yaw", 180., "偏航角");
  params.addParam<Real>("roll", -90., "滚转角");

  params.addParam<Real>("ref_length", 1, "参考长度");
  params.addParam<Real>("ref_area", 1, "参考面积");
  params.addParam<Real>("ds", 1.490116119384766e-08, "微扰量");
  return params;
}

EulerBase::EulerBase(const std::string & name, InputParameters parameters)
{
	_mach = (parameters.get<Real>("mach"));
	_gamma = (parameters.get<Real>("gamma"));

	_attack = (parameters.get<Real>("attack"));
	_sideslip = (parameters.get<Real>("sideslip"));
	_pitch = (parameters.get<Real>("pitch"));
	_yaw = (parameters.get<Real>("yaw"));
	_roll = (parameters.get<Real>("roll"));

	_attack *= libMesh::pi/180;
	_sideslip *= libMesh::pi/180;
	_pitch *= libMesh::pi/180;
	_yaw *= libMesh::pi/180;
	_roll *= libMesh::pi/180;

	_ref_length = (parameters.get<Real>("ref_length"));
	_ref_area = (parameters.get<Real>("ref_area"));

	_ds = (parameters.get<Real>("ds"));
}

Real EulerBase::pressure(Real *uh)
{
	return (_gamma-1) * (uh[4] - 0.5*(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0]);  //
}

Real EulerBase::pressureInfity()
{
	return 1.0/_gamma/_mach/_mach;
}
Real EulerBase::enthalpy(Real *uh)
{
	Real p = pressure(uh);
	return (uh[4] + p)/uh[0];
}

Real EulerBase::temperature(Real* uh)
{
	Real p = pressure(uh);
	return _gamma*_mach*_mach*p/uh[0];
}

Real EulerBase::mach_local(Real* uh)
{
	Real vel = std::sqrt(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0];
	Real c = std::sqrt(temperature(uh))/_mach;
	return vel/c;
}

Real EulerBase::maxEigenValue(Real *uh, const Point &normal)
{
	RealVectorValue vel(uh[1]/uh[0], uh[2]/uh[0], uh[3]/uh[0]);
	Real c = std::sqrt(_gamma*pressure(uh)/uh[0]);
	return std::fabs(vel*normal)+c;
}

void EulerBase::eigenValue(Real *lam, Real *uh, const Point &normal)
{
	RealVectorValue vel(uh[1]/uh[0], uh[2]/uh[0], uh[3]/uh[0]);
	Real c = std::sqrt(_gamma*pressure(uh)/uh[0]);
	lam[0] = vel*normal-c;
	lam[1] = vel*normal;
	lam[2] = vel*normal+c;
}

void EulerBase::inviscousTerm(RealVectorValue* inviscous_term, Real* uh)
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

void EulerBase::inviscousTerm(std::vector<RealVectorValue>& inviscous_term, Real* uh)
{
	inviscousTerm(&inviscous_term[0], uh);
}

Quaterniond EulerBase::bodyFromWind()
{
	Quaterniond q_attack(AngleAxisd(_attack, Vector3d::UnitY()));
	Quaterniond q_sideslip (AngleAxisd(pi-_sideslip,  Vector3d::UnitZ()));
	return q_attack*q_sideslip;
}

Quaterniond EulerBase::earthFromBody()
{
	Quaterniond q_pitch(AngleAxisd(_pitch, Vector3d::UnitY()));
	Quaterniond q_yaw(AngleAxisd(_yaw, Vector3d::UnitZ()));
	Quaterniond q_roll(AngleAxisd(_roll, Vector3d::UnitX()));
	return q_roll*q_pitch*q_yaw;
}

Quaterniond EulerBase::earthFromWind()
{
	return earthFromBody()*bodyFromWind();
}

int EulerBase::equationIndex(const std::string &var_name)
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
