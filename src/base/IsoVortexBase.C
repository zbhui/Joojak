
#include "IsoVortexBase.h"

template<>
InputParameters validParams<IsoVortexBase>()
{
  InputParameters params = validParams<EulerBase>();

  return params;
}

IsoVortexBase::IsoVortexBase(const std::string & name, InputParameters parameters):
		EulerBase(name, parameters)
{
}

Real IsoVortexBase::value(Real t, const Point& p, int eq)
{
	switch (eq) {
	case 0:
		return density(t, p);
		break;
	case 1:
		return x_momentum(t, p);
		break;
	case 2:
		return y_momentum(t, p);
		break;
	case 3:
		return z_momentum(t, p);
	case 4:
		return total_energy(t, p);
		break;
	default:
		return 0.0;
		mooseError("不可用的分量" << eq);
		break;
	}
}

Real IsoVortexBase::density(Real t, const Point &p)
{
	Real x = p(0)-t;
	Real y = p(1)-t;
	Real z = p(2)-t;
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

Real IsoVortexBase::x_momentum(Real t, const Point &p)
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

Real IsoVortexBase::y_momentum(Real t, const Point &p)
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

Real IsoVortexBase::z_momentum(Real t, const Point &p)
{
	return 0.0;
}


Real IsoVortexBase::total_energy(Real t, const Point &p)
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
