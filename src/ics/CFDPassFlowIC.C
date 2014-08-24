
#include "CFDPassFlowIC.h"

template<>
InputParameters validParams<CFDPassFlowIC>()
{
	InputParameters params = validParams<CFDInitialCondition>();
	params += validParams<CFDBase>();
    params.addParam<Real>("velocity", 1.0, "均匀流速度");
    return params;
}

CFDPassFlowIC::CFDPassFlowIC(const std::string & name, InputParameters parameters) :
    CFDInitialCondition(name, parameters),
    CFDBase(name, parameters),
    _velocity(getParam<Real>("velocity"))
{}

Real CFDPassFlowIC::density(const Point &p)
{
	return 1.0;
}

Real CFDPassFlowIC::x_momentum(const Point &p)
{
	Real u = _velocity * cos(_attack) * cos(_slide);
	Real v = _velocity * sin(_attack) * cos(_slide);
	Real w = _velocity * sin(_slide);

	return density(p)*u;
}

Real CFDPassFlowIC::y_momentum(const Point &p)
{
	Real u = _velocity * cos(_attack) * cos(_slide);
	Real v = _velocity * sin(_attack) * cos(_slide);
	Real w = _velocity * sin(_slide);

	return density(p)*v;
}

Real CFDPassFlowIC::z_momentum(const Point &p)
{
	Real u = _velocity * cos(_attack) * cos(_slide);
	Real v = _velocity * sin(_attack) * cos(_slide);
	Real w = _velocity * sin(_slide);

	return density(p)*w;
}

Real CFDPassFlowIC::total_energy(const Point &p)
{
	Real pre = 1./_gamma/_mach/_mach;
	return pre/(_gamma-1) + 0.5*density(p)*(_velocity*_velocity);
}
