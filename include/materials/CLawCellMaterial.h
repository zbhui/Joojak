
#pragma once

#include "CLawMaterial.h"
#include "CLawCellMaterialData.h"
#include "CLawMaterialData.h"

using std::vector;
class CLawProblem;

class CLawCellMaterial : public CLawMaterial
{
public:
	CLawCellMaterial(const std::string & name, InputParameters parameters);

public:
	std::vector<VariableValue*> _uh;
	std::vector<VariableGradient*> _grad_uh;
	MaterialProperty<CLawCellMaterialData >& _cell_material_data;

	virtual void computeProperties();
};

template<>
InputParameters validParams<CLawCellMaterial>();
