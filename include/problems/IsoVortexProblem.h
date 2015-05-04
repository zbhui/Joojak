
#pragma once

#include "EulerProblem.h"

class IsoVortexProblem : public EulerProblem
{
public:
	IsoVortexProblem(const std::string & name, InputParameters params);

	virtual void boundaryCondition(Real *ur, Real *ul, Point &normal, std::string bc_type);

private:

	Real density(Real t, const Point &p);
	Real momentumX(Real t, const Point &p);
	Real momentumY(Real t, const Point &p);
	Real momentumZ(Real t, const Point &p);
	Real energyTotal(Real t, const Point &p);
};

template<>
InputParameters validParams<IsoVortexProblem>();
