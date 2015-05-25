
#pragma once

#include "CFDProblem.h"

class EulerProblem : public CFDProblem
{
public:
	struct DependValue
	{
		Real pressure;
		RealVectorValue vel;
		Real viscosity;

		Real uh[10];

		void update();
	};

public:
	EulerProblem(const std::string & name, InputParameters params);

	virtual void computeCellFlux(RealGradient *flux, Real *source, Real *uh, RealGradient *duh);
	virtual void computeFaceFlux(Real* flux, RealVectorValue* lift, Real* ul, Real* ur, RealGradient* dul, RealGradient* dur, Point& normal, Real penalty);


	virtual Real physicalViscosity(Real *uh);
	virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
	virtual void viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh);
	virtual void artificialViscous(RealVectorValue* artificial_viscous, Real* uh, RealGradient *duh);
	virtual void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
	virtual void boundaryCondition(Real *ur, Real *ul, Point &normal, std::string bc_type);


private:
	void wall(Real *ur,  Real *ul, Point &normal);
	void farField(Real *ur,  Real *ul, Point &normal);
	void symmetric(Real *ur,  Real *ul, Point &normal);

	void updateDependValue(DependValue &denpend_value,  Real *uh);
};

template<>
InputParameters validParams<EulerProblem>();
