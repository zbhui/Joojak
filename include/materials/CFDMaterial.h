
#pragma once

#include "CLawMaterial.h"
#include "CFDMaterialData.h"

using std::vector;
class CFDProblem;
class CFDUserObject;

class CFDMaterial : public CLawMaterial
{
public:
	CFDMaterial(const std::string & name, InputParameters parameters);

public:
	MaterialProperty<CFDMaterialData >& _material_data;
	virtual void computeQpProperties();

private:
	  CFDProblem &_cfd_problem;
	  const CFDUserObject & _cfd_user_object;
};

template<>
InputParameters validParams<CFDMaterial>();
