
#include "KOIC.h"

template<>
InputParameters validParams<KOIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params += validParams<KOBase>();
  return params;
}

KOIC::KOIC(const std::string & name, InputParameters parameters) :
		InitialCondition(name, parameters),
		KOBase(name, parameters),
		_velocity(1)
{
	std::string var_name = _var.name();


	if(var_name == "rho")
		_eq = 0;
	if(var_name == "momentum_x")
		_eq = 1;
	if(var_name == "momentum_y")
		_eq = 2;
	if(var_name == "momentum_z")
		_eq = 3;
	if(var_name == "rhoe")
		_eq = 4;
	if(var_name == "rhok")
		_eq = 5;
	if(var_name == "rhoo")
		_eq = 6;
}

Real KOIC::value(const Point & p)
{
switch (_eq) {
	case 0:
		return density(p);
		break;
	case 1:
		return x_momentum(p);
		break;
	case 2:
		return y_momentum(p);
		break;
	case 3:
		return z_momentum(p);
		break;
	case 4:
		return total_energy(p);
		break;
	case 5:
		return turbulence_kinetic_energy(p);
		break;
	case 6:
		return turbulence_disspation_ratio(p);
		break;
	default:
		return 0.0;
		mooseError("不可用的分量" << _eq);
		break;
	}
	return 0.;
}
Real KOIC::density(const Point &p)
{
	return 1.0;
}

Real KOIC::x_momentum(const Point &p)
{
	Vector3d vel = earthFromWind()*Vector3d::UnitX();
	return density(p)*vel(0);
}

Real KOIC::y_momentum(const Point &p)
{
	Vector3d vel = earthFromWind()*Vector3d::UnitX();
	return density(p)*vel(1);
}

Real KOIC::z_momentum(const Point &p)
{
	Vector3d vel = earthFromWind()*Vector3d::UnitX();
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

Real KOIC::total_energy(const Point &p)
{
	Real pre = 1./_gamma/_mach/_mach;
	return pre/(_gamma-1) + 0.5*density(p)*(_velocity*_velocity);
}

Real KOIC::turbulence_kinetic_energy(const Point& p)
{
	return density(p)*_tu_infty;
}

Real KOIC::turbulence_disspation_ratio(const Point& p)
{
//	Real omega = _reynolds*density(p)*_tu_infty/_r_mu;
	return density(p)*log(_omega_infty);
}
