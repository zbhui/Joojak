
#pragma once

#include "FEProblem.h"

class CLawProblem : public FEProblem
{
public:
	CLawProblem(const std::string & name, InputParameters params);

	virtual int equationIndex(const std::string &var_name);
	virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
	virtual void viscousTerm(RealVectorValue *viscous_term, Real* uh, RealGradient *duh);
	virtual void sourceTerm(RealVectorValue *source_term, Real* uh, RealGradient *duh);
	virtual void fluxRiemann(Real *flux, Real *ul, Real *ur, const Point &normal);
	virtual void boundaryCondition(Real *ur, Real *ul, Point &normal, std::string bc_type);

	virtual void computeCellFlux(RealGradient *flux, Real *uh, RealGradient *duh);
	virtual void computeFaceFlux(Real *flux, RealVectorValue *lift, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur, const Point &normal, Real penalty);
	virtual void computeBoundaryFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul, Point &normal, Real penalty, std::string bc_type);

private:
	void computeLift(RealVectorValue *lift, Real *ul, Real *ur, const Point &normal);

public:
	int _n_equations;
};

template<>
InputParameters validParams<CLawProblem>();
