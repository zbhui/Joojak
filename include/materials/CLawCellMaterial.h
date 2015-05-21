
#pragma once

#include "CLawMaterial.h"

using std::vector;
class CLawProblem;

class CLawCellMaterialData
{
public:
	RealVectorValue _flux_term[10];
	RealVectorValue _flux_jacobi_variable[10][10];
	RealTensorValue _flux_jacobi_grad_variable[10][10];
	Real _source_term[10];
	Real _source_jacobi_variable[10][10];
	RealVectorValue _source_jacobi_grad_variable[10][10];
};

class CLawCellMaterial : public CLawMaterial
{
public:
	CLawCellMaterial(const std::string & name, InputParameters parameters);

public:
	MaterialProperty<CLawCellMaterialData >& _material_data;
	virtual void computeProperties();
};

template<>
InputParameters validParams<CLawCellMaterial>();
