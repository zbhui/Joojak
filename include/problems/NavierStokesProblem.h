
#pragma once

#include "CLawProblem.h"
#include "Eigen/Geometry"
using Eigen::Quaterniond;

class NavierStokesProblem : public CLawProblem
{
public:
	NavierStokesProblem(const std::string & name, InputParameters params);

	int equationIndex(const std::string &var_name);
	Real pressure(Real *uh);
	Real pressureInfity();
	Real enthalpy(Real *uh);
	Real temperature(Real  *uh);
	Real mach_local(Real *uh);
	Real acous(Real *uh);
	Real maxEigenValue(Real *uh, const Point &normal);
	void eigenValue(Real *lam, Real *uh, const Point &normal);
	virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
	virtual void inviscousTerm(std::vector<RealVectorValue> &inviscous_term, Real *uh);
	virtual Quaterniond bodyFromWind();
	virtual Quaterniond earthFromBody();
	virtual Quaterniond earthFromWind();

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

	Real _ref_length;
	Real _ref_area;
};

template<>
InputParameters validParams<NavierStokesProblem>();
