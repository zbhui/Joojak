
#include "CFDInitialCondition.h"

template<>
InputParameters validParams<CFDInitialCondition>()
{
  InputParameters params = validParams<InitialCondition>();
  return params;
}

CFDInitialCondition::CFDInitialCondition(const std::string & name, InputParameters parameters) :
    InitialCondition(name, parameters)
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
}

Real CFDInitialCondition::value(const Point & p)
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
	default:
		return 0.0;
		mooseError("不可用的分量" << _eq);
		break;
	}
}
