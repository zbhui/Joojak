
#pragma once

#include "Moose.h"
#include "Eigen/Geometry"
using Eigen::Quaterniond;

class Attitude
{
public:
	Attitude(Real attack, Real sideslip, Real pitch, Real yaw, Real roll);

public:
	Quaterniond bodyFromWind();
	Quaterniond earthFromBody();
	Quaterniond earthFromWind();

protected:
	Real _attack;			/// 攻角
	Real _sideslip;		///侧滑角
	Real _pitch;			///俯仰角
	Real _yaw;			///偏航角
	Real _roll;			///滚转角
};
