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

class EulerCellMaterial;

template<>
InputParameters validParams<EulerCellMaterial>();

/**
 * Euler流体的材料属性
 */
class EulerCellMaterial :
public Material,
public EulerBase
{
public:
	EulerCellMaterial(const std::string & name, InputParameters parameters);

protected:
	virtual void computeQpProperties();

	int _n_equations;
	/// 积分点上的变量值
	std::vector<VariableValue*> _uh;
	MaterialProperty<std::vector<RealVectorValue> > & _invis_term;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > >& _jacobi;
	VariableValue &_disp_x;
	VariableValue &_disp_y;
	VariableValue &_disp_old_x;
	VariableValue &_disp_old_y;
	void computeQpValue(Real *uh);
	void fluxTerm(RealVectorValue *flux_term, Real *uh);
};
