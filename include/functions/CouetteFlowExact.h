
#pragma once

#include "Function.h"
#include "NSBase.h"

class CouetteFlowExact;

template<>
InputParameters validParams<CouetteFlowExact>();

class CouetteFlowExact :
public Function,
public NSBase
{
public:
  CouetteFlowExact(const std::string & name, InputParameters parameters);

  Real value(Real t, const Point & p);

protected:
};

