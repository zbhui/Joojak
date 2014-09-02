
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
  Real physicalViscosity(Real *uh);
  virtual void viscousTerm(RealVectorValue *viscous_term, Real *uh, RealGradient *duh);
  virtual void stressTerm(RealTensorValue &stress_term, Real *uh, RealGradient *duh);
protected:
  Real _prandtl;
  Real _reynolds;

  Real _epsilon;
  Real _sigma;

};

