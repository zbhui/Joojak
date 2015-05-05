
#pragma once

#include "InputParameters.h"

class NavierStokesProblem;

class CouetteFlowBase
{
public:
	CouetteFlowBase(const std::string & name, InputParameters parameters);
	Real value(Real t, const Point &p, int eq);

protected:
	Real density(Real t, const Point &p);
	Real momentumX(Real t, const Point &p);
	Real momentumY(Real t, const Point &p);
	Real momentumZ(Real t, const Point &p);
	Real energyTotal(Real t, const Point &p);
	Real temperature(Real t, const Point& p);

	NavierStokesProblem & _ns_problem;
	Real _gamma;
	Real _mach;
	Real _reynolds;
	Real _prandtl;
};

template<>
InputParameters validParams<CouetteFlowBase>();
