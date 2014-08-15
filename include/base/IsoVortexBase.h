
#pragma once

#include "CFDBase.h"

class IsoVortexBase;

template<>
InputParameters validParams<IsoVortexBase>();

/**
 * CFD base class.
 */
class IsoVortexBase : public CFDBase
{
public:
	IsoVortexBase(const std::string & name, InputParameters parameters);

protected:
	Real density(const Point &p);
	Real x_momentum(const Point &p);
	Real y_momentum(const Point &p);
	Real z_momentum(const Point &p);
	Real total_energy(const Point &p);

};

