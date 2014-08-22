
#include "EulerFaceKernel.h"

template<>
InputParameters validParams<EulerFaceKernel>()
{
	  InputParameters params = validParams<DGKernel>();
	  return params;
}
EulerFaceKernel::EulerFaceKernel(const std::string & name, InputParameters parameters):
		DGKernel(name, parameters),
		_flux(getNeighborMaterialProperty<std::vector<Real> >("flux")),
		_jacobi_variable_ee(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("face_jacobi_ee")),
		_jacobi_variable_en(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("face_jacobi_en")),
		_jacobi_variable_ne(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("face_jacobi_ne")),
		_jacobi_variable_nn(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("face_jacobi_nn"))
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
	Real flux = _flux[_qp][_eq];
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
	Real r = 0;
	switch (type)
	{
	case Moose::ElementElement:
		r = _jacobi_variable_ee[_qp][_eq][_eq]*_phi[_j][_qp]*_test[_i][_qp];
		break;

	case Moose::ElementNeighbor:
		r = _jacobi_variable_en[_qp][_eq][_eq]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		break;

	case Moose::NeighborElement:
		r = _jacobi_variable_ne[_qp][_eq][_eq]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		break;

	case Moose::NeighborNeighbor:
		r = _jacobi_variable_nn[_qp][_eq][_eq]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		break;
	}

	return r;
}

Real EulerFaceKernel::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
	Real r = 0;
	switch (type)
	{
	case Moose::ElementElement:
		r = _jacobi_variable_ee[_qp][_eq][jvar]*_phi[_j][_qp]*_test[_i][_qp];
		break;

	case Moose::ElementNeighbor:
		r = _jacobi_variable_en[_qp][_eq][jvar]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		break;

	case Moose::NeighborElement:
		r = _jacobi_variable_ne[_qp][_eq][jvar]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		break;

	case Moose::NeighborNeighbor:
		r = _jacobi_variable_nn[_qp][_eq][jvar]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		break;
	}

	return r;
}
