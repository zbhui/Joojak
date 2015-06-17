
#pragma once

#include "CLawMaterial.h"
#include "CFDMaterialData.h"

using std::vector;
class CLawProblem;

class CFDMaterial : public CLawMaterial
{
public:
	CFDMaterial(const std::string & name, InputParameters parameters);

public:
	MaterialProperty<CFDMaterialData >& _material_data;
	virtual void computeQpProperties();
};

template<>
InputParameters validParams<CFDMaterial>();
