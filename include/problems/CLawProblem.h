
#pragma once

#include "FEProblem.h"
using std::vector;

class CLawCellMaterial;
class CLawFaceMaterial;
class CLawBoundaryMaterial;

class CLawProblem :
public FEProblem
{
public:
	CLawProblem(const std::string & name, InputParameters params);

	virtual void computeCellFlux(RealGradient *flux, Real *uh, RealGradient *duh);
	virtual void computeCellFlux(RealGradient *flux, Real *source, Real *uh, RealGradient *duh);
	virtual void computeFaceFlux(Real *flux, RealVectorValue *lift, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur, Point &normal, Real penalty);
	virtual void computeBoundaryFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul, Point &normal, Real penalty, std::string bc_type);
	virtual void computeCellMaterial(CLawCellMaterial& claw_cell_material);
	virtual void computeFaceMaterial(CLawFaceMaterial& claw_face_material);
	virtual void computeBoundaryMaterial(CLawBoundaryMaterial& bnd_material);

	const vector<VariableName> & getAuxVariables() {return _aux_variables;}
protected:
	virtual void computeLift(RealVectorValue *lift, Real *ul, Real *ur, Point &normal);
	virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
	virtual void viscousTerm(RealVectorValue *viscous_term, Real* uh, RealGradient *duh);
	virtual void sourceTerm(Real *source_term, Real* uh, RealGradient *duh);
	virtual void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
	virtual void boundaryCondition(Real *ur, Real *ul, Point &normal, std::string bc_type);
	virtual Real initialCondition(int eq);

	virtual void init();


public:
	vector<VariableName> _aux_variables;
	int _n_equations;

	vector<Real> _uh;
	vector<RealGradient> _duh;
	vector<RealVectorValue> _inviscous_term;
	vector<RealVectorValue> _viscous_term;
};

template<>
InputParameters validParams<CLawProblem>();
