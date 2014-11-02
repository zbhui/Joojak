
#include "SAIC.h"

template<>
InputParameters validParams<SAIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params += validParams<SABase>();
  params.addParam<Real>("velocity", 1.0, "均匀流速度");
  return params;
}

SAIC::SAIC(const std::string & name, InputParameters parameters) :
		InitialCondition(name, parameters),
		SABase(name, parameters),
	    _velocity(getParam<Real>("velocity"))
{
	_eq = equationIndex(_var.name());
}

Real SAIC::value(const Point & p)
{
switch (_eq) {
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
		return totalEnergy(p);
		break;
	case 5:
		return eddyViscosity(p);
		break;
	default:
		return 0.0;
		mooseError("不可用的分量" << _eq);
		break;
	}
	return 0.;
}
Real SAIC::density(const Point &p)
{
	return 1.0;
}

Real SAIC::momentumX(const Point &p)
{
	Vector3d vel = _velocity*(earthFromWind()*Vector3d::UnitX());
	return density(p)*vel(0);
}

Real SAIC::momentumY(const Point &p)
{
	Vector3d vel = _velocity*(earthFromWind()*Vector3d::UnitX());
	return density(p)*vel(1);
}

Real SAIC::momentumZ(const Point &p)
{
	Vector3d vel = _velocity*(earthFromWind()*Vector3d::UnitX());
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

Real SAIC::totalEnergy(const Point &p)
{
	Real pre = 1./_gamma/_mach/_mach;
	return pre/(_gamma-1) + 0.5*density(p)*(_velocity*_velocity);
}

Real SAIC::eddyViscosity(const Point& p)
{
	return density(p)*_nu_infty;
}
