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
#include "NSBase.h"

class NSBndMaterial;

template<>
InputParameters validParams<NSBndMaterial>();

class NSBndMaterial :
public Material,
public NSBase
{
public:
	NSBndMaterial(const std::string & name, InputParameters parameters);

protected:
	MooseEnum _bc_type;
	int _n_equations;

	/// 积分点上的变量值
	std::vector<VariableValue*> _ul;
	std::vector<VariableGradient*> _grad_ul;
	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > &_flux_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealGradient> > > &_flux_jacobi_grad_variable;
	MaterialProperty<std::vector<RealVectorValue> > & _penalty;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > &_penalty_jacobi_variable;

	virtual void computeQpProperties();
	virtual void resizeQpProperty();
	virtual void computeQpLeftValue(Real *ul, RealGradient *dul);
	virtual void computeQpRightValue(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);

	void fluxTerm(Real *flux, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	void penaltyTerm(RealVectorValue *penalty, Real *ul, Real *ur);

	void wall(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);
	void farField(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);
	void symmetric(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);
};
