
#pragma once

#include "CFDProblem.h"
#include "Attitude.h"

class EulerProblem : public CFDProblem
{
public:
	EulerProblem(const std::string & name, InputParameters params);

	virtual Real physicalViscosity(Real *uh);
	virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
	virtual void viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh);
	virtual void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
	virtual void boundaryCondition(Real *ur, Real *ul, Point &normal, std::string bc_type);

private:
	void wall(Real *ur,  Real *ul, Point &normal);
	void farField(Real *ur,  Real *ul, Point &normal);
	void symmetric(Real *ur,  Real *ul, Point &normal);

};

template<>
InputParameters validParams<EulerProblem>();
