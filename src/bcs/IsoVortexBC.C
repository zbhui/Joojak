//
//#include "IsoVortexBC.h"
//
//template<>
//InputParameters validParams<IsoVortexBC>()
//{
//	InputParameters params = validParams<EulerBC>();
//
//	return params;
//}
//
//IsoVortexBC::IsoVortexBC(const std::string &name, InputParameters parameters):
//		EulerBC(name, parameters)
//{
//}
//
//void IsoVortexBC::valueAtRightFace(Real *ur)
//{
//	Real ui[5];
//	valueAtLeftFace(ui);
//    Point p = _q_point[_qp];
//    ur[0] = density(_t, p);
//    ur[1] = x_momentum(_t, p);
//    ur[2] = y_momentum(_t, p);
//    ur[3] = z_momentum(_t, p);
//    ur[4] = total_energy(_t, p);
//}
//
//Real IsoVortexBC::density(Real t, const Point &p)
//{
//	Real x = p(0);
//	Real y = p(1);
//	Real z = p(2);
//
//	Real gam = 1.4, gamm1 = gam - 1, epi = 5.0;
//	Real xb, yb, r2;
//	Real rho,T;
//
//	Real PI = libMesh::pi;
//	xb = x+5-t;
//	yb = y+5-t;
//	r2 = xb * xb + yb * yb;
//
//	T = 1.0 - gamm1 * epi * epi / ( 8 * gam * PI* PI ) * exp( 1 - r2 );
//	rho = pow( T, 1 / gamm1 );
//	return rho;
//
//}
//
//Real IsoVortexBC::x_momentum(Real t, const Point &p)
//{
//	Real x = p(0);
//	Real y = p(1);
//	Real z = p(2);
//
//	Real gam=1.4, gamm1=gam-1, epi=5.;
//	Real xb, yb, r2;
//	Real rho, u, T;
//
//	Real PI = libMesh::pi;
//	xb = x+5-t;
//	yb = y+5-t;
//	r2=xb*xb+yb*yb;
//	u = 1.+epi/( 2.*PI ) * exp( 0.5 * ( 1.-r2 ) ) * (-yb );
//	T = 1. - gamm1*epi*epi/( 8*gam*PI*PI ) * exp( 1.-r2 );
//	rho = pow( T, 1/gamm1 );
//
//	return rho*u;
//}
//
//Real IsoVortexBC::y_momentum(Real t, const Point &p)
//{
//	Real x = p(0);
//	Real y = p(1);
//	Real z = p(2);
//
//	Real gam=1.4, gamm1=gam-1, epi=5.;
//	Real xb, yb, r2;
//	Real rho, v,T;
//
//	Real PI = libMesh::pi;
//	xb = x+5-t;
//	yb = y+5-t;
//	r2=xb*xb+yb*yb;
//	v = 1.+epi/( 2.*PI ) * exp( 0.5 * ( 1.-r2 ) ) * xb;
//	T = 1. - gamm1*epi*epi/( 8*gam*PI*PI ) * exp( 1.-r2 );
//	rho = pow( T, 1/gamm1 );
//
//	return rho*v;
//}
//
//Real IsoVortexBC::z_momentum(Real t, const Point &p)
//{
//	return 0.0;
//}
//
//Real IsoVortexBC::total_energy(Real t, const Point &p)
//{
//	Real x = p(0);
//	Real y = p(1);
//	Real z = p(2);
//
//	Real gam=1.4, gamm1=gam-1, epi=5.;
//	Real xb, yb, r2;
//	Real rho, u, v,T, pre;
//
//	Real PI = libMesh::pi;
//	xb = x+5-t;
//	yb = y+5-t;
//	r2=xb*xb+yb*yb;
//	u = 1.+epi/( 2.*PI ) * exp( 0.5 * ( 1.-r2 ) ) * (-yb );
//	v = 1.+epi/( 2.*PI ) * exp( 0.5 * ( 1.-r2 ) ) * xb;
//	T = 1. - gamm1*epi*epi/( 8*gam*PI*PI ) * exp( 1.-r2 );
//	rho = pow( T, 1/gamm1 );
//	pre=pow( rho, gam );
//
//	return pre/gamm1+0.5*rho * ( u*u+v*v );
//}
//
