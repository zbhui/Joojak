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

class EulerBndMaterial;

template<>
InputParameters validParams<EulerBndMaterial>();

class EulerBndMaterial :
public Material,
public EulerBase
{
public:
	EulerBndMaterial(const std::string & name, InputParameters parameters);

protected:
	MooseEnum _bc_type;
	int _n_equations;

	/// 积分点上的变量值
	std::vector<VariableValue*> _ul;
	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > &_jacobi_variable;

	virtual void computeQpProperties();
	virtual void resizeQpProperty();
	virtual void computeQpLeftValue(Real *ul);
	virtual void computeQpRightValue(Real *ur);

	void fluxRiemann(Real *flux, Real *ul, Real *ur);
	void wall(Real *ur);
	void farField(Real *ur);
	void symmetric(Real *ur);
};
