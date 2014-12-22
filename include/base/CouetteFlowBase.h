
#pragma once

#include "NSBase.h"

class CouetteFlowBase;

template<>
InputParameters validParams<CouetteFlowBase>();

class CouetteFlowBase : public NSBase
{
public:
	CouetteFlowBase(const std::string & name, InputParameters parameters);

	Real value(Real t, const Point &p, int eq);
protected:
	Real density(Real t, const Point &p);
	Real x_momentum(Real t, const Point &p);
	Real y_momentum(Real t, const Point &p);
	Real z_momentum(Real t, const Point &p);
	Real total_energy(Real t, const Point &p);
	Real temperature(const Point& p);
};

