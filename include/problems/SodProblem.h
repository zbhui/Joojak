
#pragma once

#include "EulerProblem.h"

class SodProblem : public EulerProblem
{
public:
	virtual Real initialCondition(const Point & point, int eq);
	SodProblem(const std::string & name, InputParameters params);

	virtual void computeBoundaryFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul, CLawBoundaryMaterial &bnd);

	Real valueExact(Real t, const Point &p, int eq);
private:
	Real density(Real t, const Point &p);
	Real momentumX(Real t, const Point &p);
	Real momentumY(Real t, const Point &p);
	Real momentumZ(Real t, const Point &p);
	Real energyTotal(Real t, const Point &p);
	Real pressure(Real t, const Point &p);
};

template<>
InputParameters validParams<SodProblem>();
