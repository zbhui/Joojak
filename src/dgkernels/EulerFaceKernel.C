
#include "EulerFaceKernel.h"

template<>
InputParameters validParams<EulerFaceKernel>()
{
	  InputParameters params = validParams<DGKernel>();
	  return params;
}
EulerFaceKernel::EulerFaceKernel(const std::string & name, InputParameters parameters):
		DGKernel(name, parameters),
		_invis_term(getMaterialProperty<RealVectorValue*>("inviscous")),
		_invis_term_neighbor(getNeighborMaterialProperty<RealVectorValue*>("inviscous"))
{
	std::string var_name = _var.name();

	if(var_name == "rho")
		_eq = 0;
	if(var_name == "momentum_x")
		_eq = 1;
	if(var_name == "momentum_y")
		_eq = 2;
	if(var_name == "momentum_z")
		_eq = 3;
	if(var_name == "rhoe")
		_eq = 4;
}

Real EulerFaceKernel::computeQpResidual(Moose::DGResidualType type)
{
	Real flux = 0.5*(_invis_term[_qp][_eq] + _invis_term_neighbor[_qp][_eq])*_normals[_qp];
	flux += 1*(_u[_qp]-_u_neighbor[_qp]);
	if(type == Moose::Element)
	{
		return flux * _test[_i][_qp];
	}
	if(type == Moose::Neighbor)
	{
		return -flux* _test_neighbor[_i][_qp];
	}
	mooseError("face flux error.");
	return 0.;
}

Real EulerFaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
	return 0.;
}
