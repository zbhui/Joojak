
#include "CFDBC.h"

template<>
InputParameters validParams<CFDBC>()
{
	MooseEnum bc_types("wall, far_field, symmetric, pressure_out, none", "none");  // 边界条件的类型，可以增加

	InputParameters params = validParams<IntegratedBC>();
	params += validParams<CFDBase>();

	params.addRequiredParam<MooseEnum>("bc_type", bc_types, "边界条件");
	return params;
}

CFDBC::CFDBC(const std::string & name, InputParameters parameters):
		IntegratedBC(name, parameters),
		CFDBase(name, parameters),
		_bc_type(getParam<MooseEnum>("bc_type"))
{
}




