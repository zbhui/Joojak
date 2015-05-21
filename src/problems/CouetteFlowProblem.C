
#include "CouetteFlowProblem.h"
#include "CLawBoundaryMaterial.h"

template<>
InputParameters validParams<CouetteFlowProblem>()
{
  InputParameters params = validParams<NavierStokesProblem>();

  return params;
}

CouetteFlowProblem::CouetteFlowProblem(const std::string & name, InputParameters params) :
	NavierStokesProblem(name, params)
{
}

Real CouetteFlowProblem::density(Real t, const Point &p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);

	Real tem = temperature(t, p);
	Real pre = 1./(_gamma * _mach * _mach);
	return pre * _gamma * _mach * _mach/tem;

}

Real CouetteFlowProblem::momentumX(Real t, const Point &p)
{
	Real rho = density(t, p);
	Real u = p(1)/2.0;
	return rho * u;
}

Real CouetteFlowProblem::momentumY(Real t, const Point &p)
{
	Real rho = density(t, p);
	Real v = 0.;
	return rho * v;
}

Real CouetteFlowProblem::momentumZ(Real t, const Point &p)
{
	return 0.0;
}

Real CouetteFlowProblem::energyTotal(Real t, const Point &p)
{
	Real rho = density(t, p);
	Real tem = temperature(t, p);
	Real pre = 1./(_gamma * _mach * _mach);
	Real u = p(1)/2.0;
	Real v = 0;
	Real w = 0;
	return pre/(_gamma-1)+0.5*rho * (u*u + v*v + w*w);
}

Real CouetteFlowProblem::temperature(Real t, const Point& p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);
	Real T2=0.85,T1=0.8;
	return T1 + ( T2 - T1 ) * y / 2 + 0.5 * _prandtl * (_gamma - 1) * _mach * _mach * y / 2 * ( 1 - y / 2 );
}


Real CouetteFlowProblem::valueExact(Real t, const Point& p, int eq)
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

void CouetteFlowProblem::computeBoundaryFlux(Real* flux, RealVectorValue* lift, Real* ul, RealGradient* dul, CLawBoundaryMaterial& bnd)
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

Real CouetteFlowProblem::initialCondition(const Point& point, int eq)
{
	return valueExact(0, point, eq);
}
