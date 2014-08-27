
#include "CouetteFlowExact.h"

template<>
InputParameters validParams<CouetteFlowExact>()
{
  InputParameters params = validParams<Function>();
  params += validParams<CouetteFlowBase>();
  return params;
}

CouetteFlowExact::CouetteFlowExact(const std::string & name, InputParameters parameters) :
    Function(name, parameters),
    CouetteFlowBase(name, parameters)
{}

Real
CouetteFlowExact::value(Real t, const Point & p)
{

	return  density(t, p);
}
