
#include "CFDMaterial.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<CFDMaterial>()
{
  InputParameters params = validParams<CLawMaterial>();
  return params;
}

CFDMaterial::CFDMaterial(const std::string & name, InputParameters parameter):
		CLawMaterial(name, parameter),
		_material_data(declareProperty<CFDMaterialData>("cfd_material_data")),
		_cfd_problem(static_cast<CFDProblem&>(_fe_problem))
{
}

void CFDMaterial::computeQpProperties()
{
	_material_data[_qp].reinit(_cfd_problem);
}

