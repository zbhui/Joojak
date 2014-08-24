
#pragma once

#include "CFDInitialCondition.h"
#include "CFDBase.h"

class CFDPassFlowIC;

template<>
InputParameters validParams<CFDPassFlowIC>();

class CFDPassFlowIC :
public CFDInitialCondition,
public CFDBase
{
public:

  CFDPassFlowIC(const std::string & name, InputParameters parameters);

private:
  Real _velocity;
protected:
  virtual Real density(const Point &p);
  virtual Real x_momentum(const Point &p);
  virtual Real y_momentum(const Point &p);
  virtual Real z_momentum(const Point &p);
  virtual Real total_energy(const Point &p);
};
