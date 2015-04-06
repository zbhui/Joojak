
#pragma once

#include "CFDInitialCondition.h"
//#include "NavierStokesProblem.h"

class NavierStokesProblem;

class CFDPassFlowIC : public CFDInitialCondition
{
public:

  CFDPassFlowIC(const std::string & name, InputParameters parameters);

private:
  Real _velocity;
  NavierStokesProblem &_na_problem;

protected:
  virtual Real density(const Point &p);
  virtual Real x_momentum(const Point &p);
  virtual Real y_momentum(const Point &p);
  virtual Real z_momentum(const Point &p);
  virtual Real total_energy(const Point &p);
};

template<>
InputParameters validParams<CFDPassFlowIC>();
