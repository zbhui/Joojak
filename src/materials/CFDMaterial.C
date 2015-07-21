
#include "CFDMaterial.h"
#include "CFDProblem.h"
#include "CFDUserObject.h"

template<>
InputParameters validParams<CFDMaterial>()
{
  InputParameters params = validParams<CLawMaterial>();
  params.addRequiredParam<UserObjectName>("user_object", "The name of user data object to use.");

  return params;
}

CFDMaterial::CFDMaterial(const std::string & name, InputParameters parameter):
		CLawMaterial(name, parameter),
		_material_data(declareProperty<CFDMaterialData>("cfd_material_data")),
		_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
		_cfd_user_object(getUserObject<CFDUserObject>("user_object"))
{
}

void CFDMaterial::computeQpProperties()
{
	_material_data[_qp].reinit(_cfd_problem);

	if(_bnd)
	{
		if(_neighbor)
		{

		}
	}

	CFDMaterialData cfd_data;
	cfd_data.reinit(_cfd_problem);

	cfd_data.invis_flux;

	std::vector<VariableValue*> value = _cfd_user_object.value();
	for(int eq = 0; eq < 5; ++eq)
	{
		_material_data[_qp].uh[eq] = (*_cfd_user_object.value(eq))[_qp];
	}
	_material_data[_qp].reinit(_cfd_problem);
	std::cout << _material_data[_qp].s << std::endl;
}

