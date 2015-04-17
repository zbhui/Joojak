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
#include "CLawInterface.h"

class CLawBoundaryMaterial :
public Material,
public CLawInterface
{
public:
	CLawBoundaryMaterial(const std::string & name, InputParameters parameters);

protected:
//	MooseEnum _bc_type;
	std::string _bc_type;
	int _n_equations;

	std::vector<VariableValue*> _ul;
	std::vector<VariableGradient*> _grad_ul;
	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > &_flux_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealGradient> > > &_flux_jacobi_grad_variable;
	MaterialProperty<std::vector<RealVectorValue> > & _lift;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > &_lift_jacobi_variable;

	virtual void computeQpProperties();
	virtual void resizeQpProperty();
	virtual void computeQpLeftValue(Real *ul, RealGradient *dul);
	virtual void computeQpRightValue(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul);

	void fluxTerm(Real *flux, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	void liftOperator(RealVectorValue *lift, Real *ul, Real *ur);
};

template<>
InputParameters validParams<CLawBoundaryMaterial>();
