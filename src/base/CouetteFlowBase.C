
#include "CouetteFlowBase.h"

template<>
InputParameters validParams<CouetteFlowBase>()
{
  InputParameters params = validParams<NSBase>();

  return params;
}

CouetteFlowBase::CouetteFlowBase(const std::string & name, InputParameters parameters):
		NSBase(name, parameters)
{
}

Real CouetteFlowBase::value(Real t, const Point& p, int eq)
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

Real CouetteFlowBase::density(Real t, const Point &p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);

	Real tem = temperature(p);
	Real pre = 1./(_gamma * _mach * _mach);
	return pre * _gamma * _mach * _mach/tem;

}

Real CouetteFlowBase::x_momentum(Real t, const Point &p)
{
	Real rho = density(t, p);
	Real u = p(1)/2.0;
	return rho * u;
}

Real CouetteFlowBase::y_momentum(Real t, const Point &p)
{
	Real rho = density(t, p);
	Real v = 0.;
	return rho * v;
}

Real CouetteFlowBase::z_momentum(Real t, const Point &p)
{
	return 0.0;
}

Real CouetteFlowBase::total_energy(Real t, const Point &p)
{
	Real rho = density(t, p);
	Real tem = temperature(p);
	Real pre = 1./(_gamma * _mach * _mach);
	Real u = p(1)/2.0;
	Real v = 0;
	Real w = 0;
	return pre/(_gamma-1)+0.5*rho * (u*u + v*v + w*w);
}

Real CouetteFlowBase::temperature(const Point& p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);
	Real T2=0.85,T1=0.8;
	return T1 + ( T2 - T1 ) * y / 2 + 0.5 * _prandtl * (_gamma - 1) * _mach * _mach * y / 2 * ( 1 - y / 2 );
}

