
#pragma once

#include "InitialCondition.h"
#include "KOBase.h"

class KOIC;

template<>
InputParameters validParams<KOIC>();

/**
 * 为CFD提供初始条件接口
 */
class KOIC :
public InitialCondition,
public KOBase
{
public:
	KOIC(const std::string & name, InputParameters parameters);

	virtual Real value(const Point & p);

protected:
	virtual Real density(const Point &p);
	virtual Real x_momentum(const Point &p);
	virtual Real y_momentum(const Point &p);
	virtual Real z_momentum(const Point &p);
	virtual Real total_energy(const Point &p);
	virtual Real turbulence_kinetic_energy(const Point &p);
	virtual Real turbulence_disspation_ratio(const Point &p);

	unsigned int _eq;
	Real _velocity;
};
