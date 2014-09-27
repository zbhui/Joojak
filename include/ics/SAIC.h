
#pragma once

#include "InitialCondition.h"
#include "SABase.h"

class SAIC;

template<>
InputParameters validParams<SAIC>();

class SAIC :
public InitialCondition,
public SABase
{
public:
	SAIC(const std::string & name, InputParameters parameters);

	virtual Real value(const Point & p);

protected:
	virtual Real density(const Point &p);
	virtual Real momentumX(const Point &p);
	virtual Real momentumY(const Point &p);
	virtual Real momentumZ(const Point &p);
	virtual Real totalEnergy(const Point &p);
	virtual Real eddyViscosity(const Point &p);

	unsigned int _eq;
	Real _velocity;
};
