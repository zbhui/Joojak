
#pragma once

#include "CLawProblem.h"

class NavierStokesProblem : public CLawProblem
{
public:
	NavierStokesProblem(const std::string & name, InputParameters params);

protected:
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
