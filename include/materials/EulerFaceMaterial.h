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
#include "EulerBase.h"

class EulerFaceMaterial;

template<>
InputParameters validParams<EulerFaceMaterial>();

/**
 * Euler流体的材料属性
 */
class EulerFaceMaterial :
public Material,
public EulerBase
{
public:
	EulerFaceMaterial(const std::string & name, InputParameters parameters);

protected:
	virtual void computeQpProperties();

	int _n_equations;
	/// 积分点上的变量值
	std::vector<VariableValue*> _ul;
	std::vector<VariableValue*> _ur;
	std::vector<VariableGradient*> _grad_ul;
	std::vector<VariableGradient*> _grad_ur;

	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > & _jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<Real> > > & _jacobi_variable_en;
	MaterialProperty<std::vector<std::vector<Real> > > & _jacobi_variable_ne;
	MaterialProperty<std::vector<std::vector<Real> > > & _jacobi_variable_nn;

	virtual void computeQpLeftValue(Real *ul);
	virtual void computeQpRightValue(Real *ur);
	virtual void computeQpLeftGradValue(RealGradient *ul);
	virtual void computeQpRightGradValue(RealGradient *ul);
	virtual void artificalViscousTerm(RealVectorValue *viscous_term, Real *uh, RealGradient *duh);
	Real artificalViscous(Real* uh);

	void fluxTerm(Real *flux, Real* ul, Real* ur, RealGradient *dul, RealGradient *dur);
};
