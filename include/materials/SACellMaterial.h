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

class SACellMaterial;

template<>
InputParameters validParams<SACellMaterial>();

class SACellMaterial :
public Material,
public SABase
{
public:
	SACellMaterial(const std::string & name, InputParameters parameters);

protected:
	virtual void resizeQpProperty();
	virtual void computeQpProperties();

	int _n_equations;

	VariableValue &_distance;
	/// 积分点上的变量值
	std::vector<VariableValue*> _uh;
	std::vector<VariableGradient*> _grad_uh;
	MaterialProperty<std::vector<RealVectorValue> > & _flux_term;
	MaterialProperty<std::vector<Real> > & _source_term;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > >& _flux_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealTensorValue> > >& _flux_jacobi_grad_variable;
	MaterialProperty<std::vector<std::vector<Real> > >& _source_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > >& _source_jacobi_grad_variable;

    /// 势函数
    VariableValue& _psi;
    /// 势函数梯度
    VariableGradient& _grad_psi;

	void computeQpValue(Real *uh, RealGradient *duh);
	void fluxTerm(RealVectorValue *flux_term, Real *source_term, Real *uh, RealGradient *duh);
	Real distance();
};
