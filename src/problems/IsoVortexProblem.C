
#include "IsoVortexProblem.h"
#include "CLawBoundaryMaterial.h"

template<>
InputParameters validParams<IsoVortexProblem>()
{
  InputParameters params = validParams<EulerProblem>();

  return params;
}

IsoVortexProblem::IsoVortexProblem(const std::string & name, InputParameters params) :
	EulerProblem(name, params)
{
}

Real IsoVortexProblem::density(Real t, const Point &p)
{
	Real x = p(0)-t;
	Real y = p(1)-t;
	Real z = p(2)-t;
//
	Real gam = 1.4, gamm1 = gam - 1, epi = 5.0;
	Real xb, yb, r2;
	Real rho,T;

	Real PI = libMesh::pi;
	xb = x+5;
	yb = y+5;
	r2 = xb * xb + yb * yb;

	T = 1.0 - gamm1 * epi * epi / ( 8 * gam * PI* PI ) * exp( 1 - r2 );
	rho = pow( T, 1 / gamm1 );
	return rho;
}

Real IsoVortexProblem::momentumX(Real t, const Point &p)
{
	Real x = p(0)-t;
	Real y = p(1)-t;
	Real z = p(2)-t;

	Real gam=1.4, gamm1=gam-1, epi=5.;
	Real xb, yb, r2;
	Real rho, u, T;

	Real PI = libMesh::pi;
	xb = x+5;
	yb = y+5;
	r2=xb*xb+yb*yb;
	u = 1.+epi/( 2.*PI ) * exp( 0.5 * ( 1.-r2 ) ) * (-yb );
	T = 1. - gamm1*epi*epi/( 8*gam*PI*PI ) * exp( 1.-r2 );
	rho = pow( T, 1/gamm1 );

	return rho*u;
}

Real IsoVortexProblem::momentumY(Real t, const Point &p)
{
	Real x = p(0)-t;
	Real y = p(1)-t;
	Real z = p(2)-t;

	Real gam=1.4, gamm1=gam-1, epi=5.;
	Real xb, yb, r2;
	Real rho, v,T;

	Real PI = libMesh::pi;
	xb = x+5;
	yb = y+5;
	r2=xb*xb+yb*yb;
	v = 1.+epi/( 2.*PI ) * exp( 0.5 * ( 1.-r2 ) ) * xb;
	T = 1. - gamm1*epi*epi/( 8*gam*PI*PI ) * exp( 1.-r2 );
	rho = pow( T, 1/gamm1 );

	return rho*v;
}

Real IsoVortexProblem::momentumZ(Real t, const Point &p)
{
	return 0.0;
}

Real IsoVortexProblem::valueExact(Real t, const Point& p, int eq)
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

Real IsoVortexProblem::energyTotal(Real t, const Point &p)
{
	Real x = p(0)-t;
	Real y = p(1)-t;
	Real z = p(2)-t;

	Real gam=1.4, gamm1=gam-1, epi=5.;
	Real xb, yb, r2;
	Real rho, u, v,T, pre;

	Real PI = libMesh::pi;
	xb = x+5;
	yb = y+5;
	r2=xb*xb+yb*yb;
	u = 1.+epi/( 2.*PI ) * exp( 0.5 * ( 1.-r2 ) ) * (-yb );
	v = 1.+epi/( 2.*PI ) * exp( 0.5 * ( 1.-r2 ) ) * xb;
	T = 1. - gamm1*epi*epi/( 8*gam*PI*PI ) * exp( 1.-r2 );
	rho = pow( T, 1/gamm1 );
	pre=pow( rho, gam );

	return pre/gamm1+0.5*rho * ( u*u+v*v );
}

void IsoVortexProblem::computeBoundaryFlux(Real* flux, RealVectorValue* lift, Real* ul, RealGradient* dul, CLawBoundaryMaterial& bnd)
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

Real IsoVortexProblem::initialCondition(const Point& point, int eq)
{
	return valueExact(0, point, eq);
}
