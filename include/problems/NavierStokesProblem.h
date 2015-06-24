
#pragma once

#include "CFDProblem.h"

class NavierStokesProblem : public CFDProblem
{
public:
	NavierStokesProblem(const std::string & name, InputParameters params);

	virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
	virtual void viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh);
	virtual void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
	virtual void boundaryCondition(Real *ur, Real *ul, Point &normal, std::string bc_type);
	virtual void computeBoundaryFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul, Point &normal, Real penalty, std::string bc_type);

protected:
	void isothermalWall(Real *ur,  Real *ul, Point &normal);
	void adiabaticWall(Real *ur,  Real *ul, Point &normal);
	void farField(Real *ur,  Real *ul, Point &normal);
	void farFieldHarteamann(Real *ur,  Real *ul, Point &normal);
	void symmetric(Real *ur,  Real *ul, Point &normal);

	void viscousTermAdiabatic(RealVectorValue* viscous_term, Real* uh, RealGradient *duh);

	struct DependValue
	{
		Real pressure;
		RealVectorValue vel;
		Real viscosity;

		Real uh[10];

		void update();
	};

	DependValue _value;
	DependValue _value_neighbor;
};

template<>
InputParameters validParams<NavierStokesProblem>();
