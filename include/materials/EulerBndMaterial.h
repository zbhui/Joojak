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
#include "CFDBase.h"

class EulerBndMaterial;

template<>
InputParameters validParams<EulerBndMaterial>();

class EulerBndMaterial :
public Material,
public CFDBase
{
public:
	EulerBndMaterial(const std::string & name, InputParameters parameters);

protected:
	virtual void computeQpProperties();

	int _n_equations;
	/// 积分点上的变量值
	std::vector<VariableValue*> _ul;
	std::vector<VariableValue*> _ur;

	MaterialProperty<std::vector<RealVectorValue> > & _invis_term;
	MaterialProperty<std::vector<RealVectorValue> > & _invis_term_neighbor;
	MaterialProperty<Real> & _flux_diff;

	virtual void computeQpLeftValue(Real *ul);
	virtual void computeQpRightValue(Real *ur);
};
