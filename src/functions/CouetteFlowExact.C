
#include "CouetteFlowExact.h"

template<>
InputParameters validParams<CouetteFlowExact>()
{
  InputParameters params = validParams<Function>();
  params += validParams<NSBase>();
  return params;
}

CouetteFlowExact::CouetteFlowExact(const std::string & name, InputParameters parameters) :
    Function(name, parameters),
    NSBase(name, parameters)
{}

Real
CouetteFlowExact::value(Real t, const Point & p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);

	Real T2=0.85,T1=0.8;
	Real tem = T1 + ( T2 - T1 ) * y / 2 + 0.5 * _prandtl * (_gamma - 1) * _mach * _mach * y / 2 * ( 1 - y / 2 );
	Real pre = 1./(_gamma * _mach * _mach);
	return pre * _gamma * _mach * _mach/tem;
}
