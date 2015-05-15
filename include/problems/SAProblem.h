
#pragma once

#include "NavierStokesProblem.h"

class SAProblem : public NavierStokesProblem
{
public:
	SAProblem(const std::string & name, InputParameters params);

	virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
	virtual void viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh);
	virtual void sourceTerm(Real* source_term, Real* uh, RealGradient* duh);
	virtual void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
	virtual void boundaryCondition(Real *ur, Real *ul, Point &normal, std::string bc_type);
	virtual void computeBoundaryFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul, Point &normal, Real penalty, std::string bc_type);


private:
	void isothermalWall(Real *ur,  Real *ul, Point &normal);
	void adiabaticWall(Real *ur,  Real *ul, Point &normal);
	void farField(Real *ur,  Real *ul, Point &normal);
	void symmetric(Real *ur,  Real *ul, Point &normal);

	void viscousTermAdiabatic(RealVectorValue* viscous_term, Real* uh, RealGradient *duh);

	Real _cb1, _cb2, _sigma_sa, _kappa;
	Real _cw1, _cw2, _cw3, _cv1, _cv2, _cv3;
	Real _ct1, _ct2, _ct3, _ct4;
	Real _prandtl_turb;

	Real _cw3_pow6;

	Real _nu_infty;
//	MooseEnum _bc_types;
};

template<>
InputParameters validParams<SAProblem>();
