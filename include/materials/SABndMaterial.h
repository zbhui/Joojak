/** ***********************************************************
*  @file
*  @brief
*  @author	刘  伟
*
*  Program:   Joojak
*  Copyright (c) 刘伟，张来平，2014，空气动力学国家重点实验室(SKLA)
*  All rights reserved.
*  ************************************************************
**/

#pragma once

#include "Material.h"
#include "SABase.h"

class SABndMaterial;

template<>
InputParameters validParams<SABndMaterial>();

class SABndMaterial :
public Material,
public SABase
{
public:
	SABndMaterial(const std::string & name, InputParameters parameters);

protected:
	const Real & _current_elem_volume;
	const Real & _current_side_volume;
	MooseEnum _bc_type;
	int _n_equations;

	/// 积分点上的变量值
	std::vector<VariableValue*> _ul;
	std::vector<VariableGradient*> _grad_ul;
	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > &_flux_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealGradient> > > &_flux_jacobi_grad_variable;
	MaterialProperty<std::vector<RealVectorValue> > & _penalty;
	MaterialProperty<std::vector<RealVectorValue> > & _penalty_neighbor;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > &_penalty_jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > &_penalty_jacobi_variable_ne;

	virtual void computeQpProperties();
	virtual void resizeQpProperty();
	virtual void computeQpLeftValue(Real *ul, RealGradient *dul);
	virtual void computeQpRightValue(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);

	void fluxTerm(Real *flux, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	void penaltyTerm(RealVectorValue* penalty, RealVectorValue* penalty_neighbor, Real* ul, Real* ur);

	void isothermalWall(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);
	void adiabaticWall(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);
	void farField(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);
	void farFieldRiemann(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);
	void symmetric(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);
	void pressureOut(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);


	virtual void viscousTerm(RealVectorValue *viscous_term, Real *uh, RealGradient *duh);
};
