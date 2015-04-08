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

class NSFaceMaterial;

template<>
InputParameters validParams<NSFaceMaterial>();

class NSFaceMaterial :
public Material,
public CLawInterface
{
public:
	NSFaceMaterial(const std::string & name, InputParameters parameters);

protected:
	virtual void resizeQpProperty();
	virtual void computeQpProperties();

	const Real & _current_elem_volume;
	const Real & _neighbor_elem_volume;
	const Real & _current_side_volume;
	Real _ds;
//	const int & _var_order;

	int _n_equations;
	/// 积分点上的变量值
	std::vector<VariableValue*> _ul;
	std::vector<VariableValue*> _ur;
	std::vector<VariableGradient*> _grad_ul;
	std::vector<VariableGradient*> _grad_ur;

	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_en;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_ne;
	MaterialProperty<std::vector<std::vector<Real> > > & _flux_jacobi_variable_nn;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_ee;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_en;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_ne;
	MaterialProperty<std::vector<std::vector<RealGradient> > > & _flux_jacobi_grad_variable_nn;

	MaterialProperty<std::vector<RealVectorValue> > & _penalty;
	MaterialProperty<std::vector<RealVectorValue> > & _penalty_neighbor;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _penalty_jacobi_variable_ee;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _penalty_jacobi_variable_en;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _penalty_jacobi_variable_ne;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > & _penalty_jacobi_variable_nn;

	virtual void computeQpLeftValue(Real *ul);
	virtual void computeQpRightValue(Real *ur);
	virtual void computeQpLeftGradValue(RealGradient *dul);
	virtual void computeQpRightGradValue(RealGradient *dur);

	void fluxTerm(Real *flux, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	void penaltyTerm(RealVectorValue *penalty, RealVectorValue *penalty_neighbor, Real *ul, Real *ur);
};
