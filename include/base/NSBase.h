
#pragma once

#include "EulerBase.h"

class NSBase;

template<>
InputParameters validParams<NSBase>();


class NSBase:
public EulerBase
{
public:
  NSBase(const std::string & name, InputParameters parameters);
  virtual ~NSBase(){};

protected:
  int equationIndex(const std::string &var_name);
  Real physicalViscosity(Real *uh);
  virtual void viscousTerm(RealVectorValue *viscous_term, Real *uh, RealGradient *duh);
  virtual void stressTerm(RealTensorValue &stress_term, Real *uh, RealGradient *duh);
  virtual void stressTerm(Matrix3d &stress_term, Real *uh, RealGradient *duh);
protected:
  Real _prandtl;
  Real _reynolds;

  Real _epsilon;
  Real _sigma;

};

