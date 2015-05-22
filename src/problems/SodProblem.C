
#include "SodProblem.h"
#include "CLawBoundaryMaterial.h"

template<>
InputParameters validParams<SodProblem>()
{
  InputParameters params = validParams<EulerProblem>();

  return params;
}

SodProblem::SodProblem(const std::string & name, InputParameters params) :
	EulerProblem(name, params)
{
}

Real SodProblem::density(Real t, const Point &p)
{
	Real x = p(0);
    if (x < 0.5)
		return 1;
    else
    	return 0.125;
}

Real SodProblem::momentumX(Real t, const Point &p)
{
	Real x = p(0);
	Real u = 0;
	Real rho = density(t, p);
    if (x < 0.5)
		u = 0;
    else
    	u = 0;

	return rho*u;
}

Real SodProblem::momentumY(Real t, const Point &p)
{
	return 0.0;
}

Real SodProblem::momentumZ(Real t, const Point &p)
{
	return 0.0;
}

Real SodProblem::energyTotal(Real t, const Point &p)
{
	RealVectorValue momentum(momentumX(t, p), momentumY(t, p), momentumZ(t, p));
	Real rho = density(t, p);
	Real pre = pressure(t, p);

	return pre/(_gamma-1) +0.5*momentum.size_sq()/rho;
}

Real SodProblem::pressure(Real t, const Point &p)
{
	Real x = p(0);
    if (x < 0.5)
		return 1;
    else
    	return 0.1;
}
void SodProblem::computeBoundaryFlux(Real* flux, RealVectorValue* lift, Real* ul, RealGradient* dul, CLawBoundaryMaterial& bnd)
{
	Point normal = bnd.normals()[_qp];
	Real penalty = bnd.penalty();
	RealGradient dur[10];
	Real ur[10];
	for (int eq = 0; eq < _n_equations; ++eq)
		dur[eq] = dul[eq];

	std::string bc_type = bnd.getBCType();
	Point q_point = bnd.qpoints()[_qp];
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = valueExact(time(), q_point, eq);

	computeFaceFlux(flux, lift, ul, ur, dul, dur, normal, penalty);
}

Real SodProblem::initialCondition(const Point& point, int eq)
{
	return valueExact(0, point, eq);
}

Real SodProblem::valueExact(Real t, const Point& p, int eq)
{
	switch (eq) {
	case 0:
		return density(t, p);
		break;
	case 1:
		return momentumX(t, p);
		break;
	case 2:
		return momentumY(t, p);
		break;
	case 3:
		return momentumZ(t, p);
	case 4:
		return energyTotal(t, p);
		break;
	default:
		return 0.0;
		mooseError("不可用的分量" << eq);
		break;
	}
}
