
#include "IsoVortexBase.h"

template<>
InputParameters validParams<IsoVortexBase>()
{
  InputParameters params = validParams<CFDBase>();

  return params;
}

IsoVortexBase::IsoVortexBase(const std::string & name, InputParameters parameters):
		CFDBase(name, parameters)
{
}

Real IsoVortexBase::density(const Point &p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);
//
//	Real pi = 3.1415926535;
//	return std::sin(2*pi*(x+y+z));

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

Real IsoVortexBase::x_momentum(const Point &p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);

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

Real IsoVortexBase::y_momentum(const Point &p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);

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

Real IsoVortexBase::z_momentum(const Point &p)
{
	return 0.0;
}

Real IsoVortexBase::total_energy(const Point &p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);

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
