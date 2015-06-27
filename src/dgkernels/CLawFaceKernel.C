
#include "CLawFaceKernel.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawFaceKernel>()
{
	  InputParameters params = validParams<DGKernel>();
	  params.addRequiredParam<int>("component", "DGKernel component");
	  return params;
}
CLawFaceKernel::CLawFaceKernel(const std::string & name, InputParameters parameters):
		DGKernel(name, parameters),
		_eq(getParam<int>("component")),
		_face(getNeighborMaterialProperty<CLawFaceMaterialData>("face_material_data"))
{
}

Real CLawFaceKernel::computeQpResidual(Moose::DGResidualType type)
{
	switch (type)
	{
	case Moose::Element:
		return  _face[_qp]._flux[_eq] * _test[_i][_qp] + _face[_qp]._lift[_eq]*_grad_test[_i][_qp];
		break;
	case Moose::Neighbor:
		return -_face[_qp]._flux[_eq] * _test_neighbor[_i][_qp] + _face[_qp]._lift[_eq]*_grad_test_neighbor[_i][_qp];
		break;
	}

	mooseError("face flux error.");
	return 0.;
}

Real CLawFaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
	return computeQpJacobian(_eq, _eq, type);
}

Real CLawFaceKernel::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
	return computeQpJacobian(_eq, jvar, type);

}

Real CLawFaceKernel::computeQpJacobian(int p, int q, Moose::DGJacobianType type)
{
	Real r = 0;
	switch (type)
	{
	case Moose::ElementElement:
		r = _face[_qp]._flux_jacobi_variable_ee[p][q]*_phi[_j][_qp]*_test[_i][_qp];
		r += _face[_qp]._flux_jacobi_grad_variable_ee[p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
		r += _face[_qp]._lift_jacobi_variable_ee[p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::ElementNeighbor:
		r = _face[_qp]._flux_jacobi_variable_en[p][q]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _face[_qp]._flux_jacobi_grad_variable_en[p][q]*_grad_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _face[_qp]._lift_jacobi_variable_en[p][q]*_grad_test[_i][_qp]*_phi_neighbor[_j][_qp];
		break;

	case Moose::NeighborElement:
		r = -_face[_qp]._flux_jacobi_variable_ee[p][q]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += -_face[_qp]._flux_jacobi_grad_variable_ee[p][q]*_grad_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += _face[_qp]._lift_jacobi_variable_ee[p][q]*_grad_test_neighbor[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::NeighborNeighbor:
		r = -_face[_qp]._flux_jacobi_variable_en[p][q]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += -_face[_qp]._flux_jacobi_grad_variable_en[p][q]*_grad_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += _face[_qp]._lift_jacobi_variable_en[p][q]*_grad_test_neighbor[_i][_qp]*_phi_neighbor[_j][_qp];
		break;
	}

	return r;
}
