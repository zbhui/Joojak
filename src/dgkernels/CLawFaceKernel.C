
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
		_flux(getNeighborMaterialProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable_ee(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_ee")),
		_flux_jacobi_variable_en(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_en")),
		_flux_jacobi_grad_variable_ee(getNeighborMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_ee")),
		_flux_jacobi_grad_variable_en(getNeighborMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_en")),

		_lift(getNeighborMaterialProperty<std::vector<RealVectorValue> >("lift")),
		_lift_jacobi_variable_ee(getNeighborMaterialProperty<std::vector<std::vector<RealVectorValue> > >("lift_jacobi_variable_ee")),
		_lift_jacobi_variable_en(getNeighborMaterialProperty<std::vector<std::vector<RealVectorValue> > >("lift_jacobi_variable_en")),
		_eq(getParam<int>("component")),
		_face_material_data(getNeighborMaterialProperty<CLawFaceMaterialData>("face_material_data"))
{
}

Real CLawFaceKernel::computeQpResidual(Moose::DGResidualType type)
{
	switch (type)
	{
	case Moose::Element:
		return  _face_material_data[_qp]._flux[_eq] * _test[_i][_qp] ;//+ _lift[_qp][_eq]*_grad_test[_i][_qp];
//		return  _flux[_qp][_eq] * _test[_i][_qp] + _lift[_qp][_eq]*_grad_test[_i][_qp];
		break;
	case Moose::Neighbor:
		return -_face_material_data[_qp]._flux[_eq] * _test_neighbor[_i][_qp] ;//+ _lift[_qp][_eq]*_grad_test_neighbor[_i][_qp];
//		return -_flux[_qp][_eq] * _test_neighbor[_i][_qp] + _lift[_qp][_eq]*_grad_test_neighbor[_i][_qp];
		break;
	}

	mooseError("face flux error.");
	return 0.;
}

Real CLawFaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
	return 0;
	return computeQpJacobian(_eq, _eq, type);
}

Real CLawFaceKernel::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
	return 0;
	return computeQpJacobian(_eq, jvar, type);

}

Real CLawFaceKernel::computeQpJacobian(int p, int q, Moose::DGJacobianType type)
{
	Real r = 0;
	switch (type)
	{
	case Moose::ElementElement:
		r = _flux_jacobi_variable_ee[_qp][p][q]*_phi[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_ee[_qp][p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
		r += _lift_jacobi_variable_ee[_qp][p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::ElementNeighbor:
		r = _flux_jacobi_variable_en[_qp][p][q]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_en[_qp][p][q]*_grad_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _lift_jacobi_variable_en[_qp][p][q]*_grad_test[_i][_qp]*_phi_neighbor[_j][_qp];
		break;

	case Moose::NeighborElement:
		r = -_flux_jacobi_variable_ee[_qp][p][q]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += -_flux_jacobi_grad_variable_ee[_qp][p][q]*_grad_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += _lift_jacobi_variable_ee[_qp][p][q]*_grad_test_neighbor[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::NeighborNeighbor:
		r = -_flux_jacobi_variable_en[_qp][p][q]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += -_flux_jacobi_grad_variable_en[_qp][p][q]*_grad_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += _lift_jacobi_variable_en[_qp][p][q]*_grad_test_neighbor[_i][_qp]*_phi_neighbor[_j][_qp];
		break;
	}

	return r;
}
