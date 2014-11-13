
#pragma once

#include "MooseVariableBase.h"
#include "MooseObject.h"
#include "Eigen/Dense"
#include "Eigen/Geometry"
using namespace Eigen;
typedef Eigen::Matrix<Real, 5, 5> Matrix5x5;

class EulerBase;

template<>
InputParameters validParams<EulerBase>();

/**
 * CFD base class.
 */
class EulerBase
{
public:
	EulerBase(const std::string & name, InputParameters parameters);
  virtual ~EulerBase(){};

protected:
  Real pressure(Real *uh);
  Real pressureInfity();
  Real enthalpy(Real *uh);
  Real temperature(Real  *uh);
  Real mach_local(Real *uh);
  Real acous(Real *uh);
  Real maxEigenValue(Real *uh, const Point &normal);
  void eigenValue(Real *lam, Real *uh, const Point &normal);
  virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
  virtual void inviscousTerm(std::vector<RealVectorValue> &inviscous_term, Real *uh);
  virtual Quaterniond bodyFromWind();
  virtual Quaterniond earthFromBody();
  virtual Quaterniond earthFromWind();
protected:

  Real _gamma;
  Real _mach;

  Real _attack;			/// 攻角
  Real _sideslip;		///侧滑角
  Real _pitch;			///俯仰角
  Real _yaw;			///偏航角
  Real _roll;			///滚转角

  Real _ref_length;
  Real _ref_area;

  Real _ds;
};

