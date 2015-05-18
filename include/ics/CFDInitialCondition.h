
#pragma once

#include "MultiInitialCondition.h"

class CFDInitialCondition : public MultiInitialCondition
{
public:
	CFDInitialCondition(const std::string & name, InputParameters parameters);

	virtual Real value(int component, const Point & p);

protected:
	virtual Real density(const Point &p) = 0;
	virtual Real momentumX(const Point &p) = 0;
	virtual Real momentumY(const Point &p) = 0;
	virtual Real momentumZ(const Point &p) = 0;
	virtual Real energyTotal(const Point &p) = 0;
	virtual Real eddyViscoisty(const Point &p) {return 0;};

};

template<>
InputParameters validParams<CFDInitialCondition>();
