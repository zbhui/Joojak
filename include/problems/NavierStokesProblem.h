
#pragma once

#include "CLawProblem.h"
#include "Eigen/Geometry"
using Eigen::Quaterniond;

class NavierStokesProblem : public CLawProblem
{
public:
	NavierStokesProblem(const std::string & name, InputParameters params);

	Real physicalViscosity(Real *uh);
	Real pressure(Real *uh);
	Real pressureInfity();
	Real enthalpy(Real *uh);
	Real temperature(Real  *uh);
	Real localMach(Real *uh);
	Real maxEigenValue(Real *uh, const Point &normal);
	void eigenValue(Real *lam, Real *uh, const Point &normal);

	virtual void inviscousTerm(std::vector<RealVectorValue> &inviscous_term, Real *uh);
	virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
	virtual void viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh);
	virtual void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
	virtual void boundaryCondition(Real *ur, Real *ul, Point &normal, std::string bc_type);
	virtual void computeBoundaryFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul, Point &normal, Real penalty, std::string bc_type);

	virtual Quaterniond bodyFromWind();
	virtual Quaterniond earthFromBody();
	virtual Quaterniond earthFromWind();

private:
	void isothermalWall(Real *ur,  Real *ul, Point &normal);
	void adiabaticWall(Real *ur,  Real *ul, Point &normal);
	void farField(Real *ur,  Real *ul, Point &normal);

	void viscousTermAdiabatic(RealVectorValue* viscous_term, Real* uh, RealGradient *duh);

public:
	Real _mach;
	Real _gamma;
	Real _reynolds;
	Real _prandtl;

	Real _attack;			/// 攻角
	Real _sideslip;		///侧滑角
	Real _pitch;			///俯仰角
	Real _yaw;			///偏航角
	Real _roll;			///滚转角

	Real _ref_length;
	Real _ref_area;

//	MooseEnum _bc_types;
};

template<>
InputParameters validParams<NavierStokesProblem>();
