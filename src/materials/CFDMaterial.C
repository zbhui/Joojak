
#include "CFDMaterial.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CFDMaterial>()
{
  InputParameters params = validParams<CLawMaterial>();
  return params;
}

CFDMaterial::CFDMaterial(const std::string & name, InputParameters parameter):
		CLawMaterial(name, parameter),
		_material_data(declareProperty<CFDMaterialData>("cfd_material_data"))
{
}

void CFDMaterial::computeQpProperties()
{
//	_claw_problem.computeCellMaterial(*this);
}

