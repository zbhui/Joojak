
#pragma once

#include "InitialCondition.h"

class CFDInitialCondition;

template<>
InputParameters validParams<CFDInitialCondition>();

/**
 * 为CFD提供初始条件接口
 */
class CFDInitialCondition :
	public InitialCondition
{
public:
	CFDInitialCondition(const std::string & name, InputParameters parameters);

	virtual Real value(const Point & p);

protected:
	virtual Real density(const Point &p) = 0;
	virtual Real x_momentum(const Point &p) = 0;
	virtual Real y_momentum(const Point &p) = 0;
	virtual Real z_momentum(const Point &p) = 0;
	virtual Real total_energy(const Point &p) = 0;

	unsigned int _eq;
};
