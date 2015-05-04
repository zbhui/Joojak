
#include "IsoVortexProblem.h"

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

void IsoVortexProblem::boundaryCondition(Real *ur, Real *ul, Point &normal, std::string bc_type)
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
