
#include "CFDProblem.h"

#include "MooseApp.h"
using namespace Eigen;

template<>
InputParameters validParams<CFDProblem>()
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
  params.addParam<Real>("init_vel", 1.0, "初始流体速度");
  return params;
}

CFDProblem::CFDProblem(const std::string & name, InputParameters params) :
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

	_attitude(Attitude(_attack, _sideslip, _pitch, _yaw, _roll)),
    _velocity(getParam<Real>("init_vel"))
//	_bc_types(MooseEnum("isothermal_wall adiabatic_wall far_field symmetric pressure_out none", "none"))
{
}

Real CFDProblem::pressure(Real *uh)
{
	return (_gamma-1) * (uh[4] - 0.5*(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0]);  //
}

Real CFDProblem::pressureInfity()
{
	return 1.0/_gamma/_mach/_mach;
}
Real CFDProblem::enthalpy(Real *uh)
{
	Real p = pressure(uh);
	return (uh[4] + p)/uh[0];
}

Real CFDProblem::temperature(Real* uh)
{
	Real p = pressure(uh);
	return _gamma*_mach*_mach*p/uh[0];
}

Real CFDProblem::localMach(Real* uh)
{
	Real vel = std::sqrt(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0];
	Real c = std::sqrt(temperature(uh))/_mach;
	return vel/c;
}

Real CFDProblem::maxEigenValue(Real *uh, const Point &normal)
{
	RealVectorValue vel(uh[1]/uh[0], uh[2]/uh[0], uh[3]/uh[0]);
	Real c = std::sqrt(_gamma*pressure(uh)/uh[0]);
	return std::fabs(vel*normal)+c;
}

void CFDProblem::eigenValue(Real *lam, Real *uh, const Point &normal)
{
	RealVectorValue vel(uh[1]/uh[0], uh[2]/uh[0], uh[3]/uh[0]);
	Real c = std::sqrt(_gamma*pressure(uh)/uh[0]);
	lam[0] = vel*normal-c;
	lam[1] = vel*normal;
	lam[2] = vel*normal+c;
}

Real CFDProblem::physicalViscosity(Real* uh)
{
	return 1.0;
}

void CFDProblem::stressTerm(RealTensorValue &tau, Real* uh, RealGradient* duh)
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

void CFDProblem::heatFluxTerm(RealVectorValue& heat_flux, Real* uh, RealGradient* duh)
{
}

Real CFDProblem::computeAuxValue(std::string var_name, Real* uh)
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

	mooseError("未知的辅助变量");
	return 0;
}

Real CFDProblem::initialCondition(const Point& p, int eq)
{
	switch (eq) {
		case 0:
			return density(p);
			break;
		case 1:
			return momentumX(p);
			break;
		case 2:
			return momentumY(p);
			break;
		case 3:
			return momentumZ(p);
			break;
		case 4:
			return energyTotal(p);
			break;
		default:
			return 0.0;
			break;
	}
}

Real CFDProblem::density(const Point &p)
{
	return 1.0;
}

Real CFDProblem::momentumX(const Point &p)
{
	Vector3d vel = _velocity*(_attitude.earthFromWind()*Vector3d::UnitX());
	return density(p)*vel(0);
}

Real CFDProblem::momentumY(const Point &p)
{
	Vector3d vel = _velocity*(_attitude.earthFromWind()*Vector3d::UnitX());
	return density(p)*vel(1);
}

Real CFDProblem::momentumZ(const Point &p)
{
	Vector3d vel = _velocity*(_attitude.earthFromWind()*Vector3d::UnitX());
	if(_mesh.dimension() == 2)
		return 0.;
	else if(_mesh.dimension() == 3)
		return density(p)*vel(2);
	else
	{
		mooseError("一维问题此处需要调试");
		return 0.;
	}
}

Real CFDProblem::energyTotal(const Point &p)
{
	Real pre = 1./_gamma/_mach/_mach;
	return pre/(_gamma-1) + 0.5*density(p)*(_velocity*_velocity);
}
