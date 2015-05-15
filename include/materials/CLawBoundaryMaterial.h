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

#include "CLawMaterial.h"

class CLawProblem;
using std::vector;

class CLawBoundaryMaterialData
{
public:
	Real _flux[10];
	Real _flux_jacobi_variable[10][10];
	RealGradient _flux_jacobi_grad_variable[10][10];
	RealVectorValue _lift[10];
	RealVectorValue _lift_jacobi_variable[10][10];
};

class CLawBoundaryMaterial : public CLawMaterial
{
public:
	CLawBoundaryMaterial(const std::string & name, InputParameters parameters);

public:
	std::string _bc_type;
	Real _ds;
	Real _sigma;
	Real _epsilon;

	std::vector<VariableValue*> _ul;
	std::vector<VariableGradient*> _grad_ul;
	MaterialProperty<CLawBoundaryMaterialData> &_material_data;

	virtual void computeProperties();
public:
	const std::string &getBCType(){return _bc_type;}
	Real penalty();

};

template<>
InputParameters validParams<CLawBoundaryMaterial>();
