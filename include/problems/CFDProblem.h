
#pragma once

#include "CLawProblem.h"
#include "Attitude.h"
#include "Eigen/Geometry"
using Eigen::Quaterniond;

class CFDProblem : public CLawProblem
{
public:
	CFDProblem(const std::string & name, InputParameters params);

	virtual Real physicalViscosity(Real *uh);
	virtual void stressTerm(RealTensorValue &tau, Real* uh, RealGradient* duh);
	virtual Real pressure(Real *uh);
	virtual Real pressureInfity();
	virtual Real enthalpy(Real *uh);
	virtual Real temperature(Real  *uh);
	virtual Real localMach(Real *uh);
	virtual Real maxEigenValue(Real *uh, const Point &normal);
	virtual void eigenValue(Real *lam, Real *uh, const Point &normal);

public:
	Real _mach;
	Real _gamma;
	Real _reynolds;
	Real _prandtl;

	Real _attack;			/// 攻角
	Real _sideslip;		///侧滑角
	Real _pitch;			///俯仰角
	Real _yaw;			///偏航角
	Real _roll;			///滚转角

	Attitude _attitude;
//	MooseEnum _bc_types;
};

template<>
InputParameters validParams<CFDProblem>();
