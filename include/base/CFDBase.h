
#pragma once

#include "MooseVariableBase.h"
#include "MooseObject.h"
#include "Eigen/Dense"
typedef Eigen::Matrix<Real, 5, 5> Matrix5x5;
typedef RealVectorValue RealVector10[10];

class CFDBase;

template<>
InputParameters validParams<CFDBase>();

/**
 * CFD base class.
 */
class CFDBase
{
public:
	CFDBase(const std::string & name, InputParameters parameters);
  virtual ~CFDBase(){};

protected:
//  virtual void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
//  virtual void fluxViscous(Real *flux, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur, Point &normal);

  Real pressure(Real *uh);
  Real enthalpy(Real *uh);
  Real temperature(Real  *uh);
  Real mach_local(Real *uh);
  Real acous(Real *uh);
  Real physicalViscosity(Real *uh);


  virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
  virtual void inviscousTerm(std::vector<RealVectorValue> &inviscous_term, Real *uh);
  virtual void viscousTerm(RealVectorValue *viscous_term, Real *uh, RealGradient *duh);

//	void liftOperator(Real *lift, Real *ul, Real *ur, Point &normal);
//	Real penaltyOperator(Real *uh, Real *delta_uh, RealGradient &grad_phi, Point &normal, int eq);
protected:

//  int _n_equation;

/// Required parameters
  Real _gamma;
  Real _prandtl;
  Real _reynolds;
  Real _mach;

  Real _attack;
  Real _slide;

  Real _epsilon;
  Real _sigma;
};

