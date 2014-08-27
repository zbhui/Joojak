
#pragma once

#include "Function.h"
#include "CouetteFlowBase.h"

class CouetteFlowExact;

template<>
InputParameters validParams<CouetteFlowExact>();

class CouetteFlowExact :
public Function,
public CouetteFlowBase
{
public:
  CouetteFlowExact(const std::string & name, InputParameters parameters);

  Real value(Real t, const Point & p);

protected:
};

