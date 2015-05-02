
#include "CFDPassFlowIC.h"
#include "CFDProblem.h"
#include "Eigen/Geometry"
using namespace Eigen;

template<>
InputParameters validParams<CFDPassFlowIC>()
{
	InputParameters params = validParams<CFDInitialCondition>();
    params.addParam<Real>("velocity", 1.0, "均匀流速度");
    return params;
}

CFDPassFlowIC::CFDPassFlowIC(const std::string & name, InputParameters parameters) :
    CFDInitialCondition(name, parameters),
    _velocity(getParam<Real>("velocity")),
	_cfd_problem(static_cast<CFDProblem&>(_fe_problem))
{
}

Real CFDPassFlowIC::density(const Point &p)
{
	return 1.0;
}

Real CFDPassFlowIC::x_momentum(const Point &p)
{
	Vector3d vel = _velocity*(_cfd_problem._attitude.earthFromWind()*Vector3d::UnitX());
	return density(p)*vel(0);
}

Real CFDPassFlowIC::y_momentum(const Point &p)
{
	Vector3d vel = _velocity*(_cfd_problem._attitude.earthFromWind()*Vector3d::UnitX());
	return density(p)*vel(1);
}

Real CFDPassFlowIC::z_momentum(const Point &p)
{
	Vector3d vel = _velocity*(_cfd_problem._attitude.earthFromWind()*Vector3d::UnitX());
	if(_current_elem->dim() == 2)
		return 0.;
	else if(_current_elem->dim() == 3)
		return density(p)*vel(2);
	else
	{
		mooseError("一维问题此处需要调试");
		return 0.;
	}
}

Real CFDPassFlowIC::total_energy(const Point &p)
{
	Real pre = 1./_cfd_problem._gamma/_cfd_problem._mach/_cfd_problem._mach;
	return pre/(_cfd_problem._gamma-1) + 0.5*density(p)*(_velocity*_velocity);
}
