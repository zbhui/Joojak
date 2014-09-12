
#include "CFDPassFlowIC.h"

template<>
InputParameters validParams<CFDPassFlowIC>()
{
	InputParameters params = validParams<CFDInitialCondition>();
	params += validParams<EulerBase>();
    params.addParam<Real>("velocity", 1.0, "均匀流速度");
    return params;
}

CFDPassFlowIC::CFDPassFlowIC(const std::string & name, InputParameters parameters) :
    CFDInitialCondition(name, parameters),
    EulerBase(name, parameters),
    _velocity(getParam<Real>("velocity"))
{}

Real CFDPassFlowIC::density(const Point &p)
{
	return 1.0;
}

Real CFDPassFlowIC::x_momentum(const Point &p)
{
	Real u = _velocity * cos(_attack) * cos(_sideslip);
	Real v = _velocity * sin(_attack) * cos(_sideslip);
	Real w = _velocity * sin(_sideslip);

	return density(p)*u;
}

Real CFDPassFlowIC::y_momentum(const Point &p)
{
	Real u = _velocity * cos(_attack) * cos(_sideslip);
	Real v = _velocity * sin(_attack) * cos(_sideslip);
	Real w = _velocity * sin(_sideslip);

	return density(p)*v;
}

Real CFDPassFlowIC::z_momentum(const Point &p)
{
	Real u = _velocity * cos(_attack) * cos(_sideslip);
	Real v = _velocity * sin(_attack) * cos(_sideslip);
	Real w = _velocity * sin(_sideslip);

	return density(p)*w;
}

Real CFDPassFlowIC::total_energy(const Point &p)
{
	Real pre = 1./_gamma/_mach/_mach;
	return pre/(_gamma-1) + 0.5*density(p)*(_velocity*_velocity);
}
