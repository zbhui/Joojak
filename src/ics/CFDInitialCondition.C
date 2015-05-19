
#include "CFDInitialCondition.h"

template<>
InputParameters validParams<CFDInitialCondition>()
{
  InputParameters params = validParams<MultiInitialCondition>();
  return params;
}

CFDInitialCondition::CFDInitialCondition(const std::string & name, InputParameters parameters) :
	MultiInitialCondition(name, parameters)
{
}

Real CFDInitialCondition::value(int component, const Point & p)
{
	std::string var_name = _var.name();
switch (component) {
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
		return energyTotal(p);
		break;
	case 5:
		return eddyViscoisty(p);
		break;
	default:
		return 0.0;

		mooseError("不可用的分量" << component );
		break;
	}
}
