
#pragma once

#include "CFDInitialCondition.h"
//#include "NavierStokesProblem.h"

class CFDProblem;

class CFDPassFlowIC : public CFDInitialCondition
{
public:

  CFDPassFlowIC(const std::string & name, InputParameters parameters);

private:
  Real _velocity;
  CFDProblem &_cfd_problem;

protected:
  virtual Real density(const Point &p);
  virtual Real momentumX(const Point &p);
  virtual Real momentumY(const Point &p);
  virtual Real momentumZ(const Point &p);
  virtual Real energyTotal(const Point &p);
};

template<>
InputParameters validParams<CFDPassFlowIC>();
