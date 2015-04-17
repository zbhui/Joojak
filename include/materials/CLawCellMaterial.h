
#pragma once

#include "Material.h"
#include "CLawInterface.h"
#include "CLawCellMaterialData.h"

class CLawCellMaterial :
public Material,
public CLawInterface
{
public:
	CLawCellMaterial(const std::string & name, InputParameters parameters);

protected:
	virtual void computeQpProperties();
	Real _ds;

	std::vector<VariableValue*> _uh;
	std::vector<VariableGradient*> _grad_uh;
	MaterialProperty<CLawCellMaterialData >& _cell_material_data;

	void fillQpValue();
	void computeMaterialData();
};

template<>
InputParameters validParams<CLawCellMaterial>();
