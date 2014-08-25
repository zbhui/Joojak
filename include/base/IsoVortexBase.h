
#pragma once

#include "EulerBase.h"

class IsoVortexBase;

template<>
InputParameters validParams<IsoVortexBase>();

/**
 * CFD base class.
 */
class IsoVortexBase : public EulerBase
{
public:
	IsoVortexBase(const std::string & name, InputParameters parameters);

	Real value(Real t, const Point &p, int eq);
protected:
	Real density(Real t, const Point &p);
	Real x_momentum(Real t, const Point &p);
	Real y_momentum(Real t, const Point &p);
	Real z_momentum(Real t, const Point &p);
	Real total_energy(Real t, const Point &p);

};

