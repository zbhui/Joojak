
#pragma once

#include "CLawInterface.h"

class EulerProblem;

class IsoVortexBase
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

	EulerProblem & _euler_problem;
	Real _gamma;
	Real _gamm1;
	Real _mach;
};

template<>
InputParameters validParams<IsoVortexBase>();
