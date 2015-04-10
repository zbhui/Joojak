
#include "CLawFaceKernel.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawFaceKernel>()
{
	  InputParameters params = validParams<DGKernel>();
	  return params;
}
CLawFaceKernel::CLawFaceKernel(const std::string & name, InputParameters parameters):
		DGKernel(name, parameters),
		CLawInterface(parameters),
		_flux(getNeighborMaterialProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable_ee(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_ee")),
		_flux_jacobi_variable_en(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_en")),
		_flux_jacobi_grad_variable_ee(getNeighborMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_ee")),
		_flux_jacobi_grad_variable_en(getNeighborMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_en")),

		_penalty(getNeighborMaterialProperty<std::vector<RealVectorValue> >("penalty")),
		_penalty_jacobi_variable_ee(getNeighborMaterialProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_ee")),
		_penalty_jacobi_variable_en(getNeighborMaterialProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_en"))
{
	_eq = _claw_problem.equationIndex(_var.name());
	_epsilon = 1;
}

Real CLawFaceKernel::computeQpResidual(Moose::DGResidualType type)
{
	switch (type)
	{
	case Moose::Element:
		return  _flux[_qp][_eq] * _test[_i][_qp] + _epsilon * _penalty[_qp][_eq]*_grad_test[_i][_qp];
		break;
	case Moose::Neighbor:
		return -_flux[_qp][_eq] * _test_neighbor[_i][_qp] + _epsilon*_penalty[_qp][_eq]*_grad_test_neighbor[_i][_qp];
		break;
	}

	mooseError("face flux error.");
	return 0.;
}

Real CLawFaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
	Real r = 0;
	int p(_eq), q(_eq);
	switch (type)
	{
	case Moose::ElementElement:
		r = _flux_jacobi_variable_ee[_qp][p][q]*_phi[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_ee[_qp][p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_ee[_qp][p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::ElementNeighbor:
		r = _flux_jacobi_variable_en[_qp][p][q]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_en[_qp][p][q]*_grad_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_en[_qp][p][q]*_grad_test[_i][_qp]*_phi_neighbor[_j][_qp];
		break;

	case Moose::NeighborElement:
		r = -_flux_jacobi_variable_ee[_qp][p][q]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += -_flux_jacobi_grad_variable_ee[_qp][p][q]*_grad_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_ee[_qp][p][q]*_grad_test_neighbor[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::NeighborNeighbor:
		r = -_flux_jacobi_variable_en[_qp][p][q]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += -_flux_jacobi_grad_variable_en[_qp][p][q]*_grad_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_en[_qp][p][q]*_grad_test_neighbor[_i][_qp]*_phi_neighbor[_j][_qp];
		break;
	}

	return r;
}

Real CLawFaceKernel::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
	Real r = 0;
	int p(_eq), q(jvar);
	switch (type)
	{
	case Moose::ElementElement:
		r = _flux_jacobi_variable_ee[_qp][p][q]*_phi[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_ee[_qp][p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_ee[_qp][p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::ElementNeighbor:
		r = _flux_jacobi_variable_en[_qp][p][q]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_en[_qp][p][q]*_grad_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_en[_qp][p][q]*_grad_test[_i][_qp]*_phi_neighbor[_j][_qp];
		break;

	case Moose::NeighborElement:
		r = -_flux_jacobi_variable_ee[_qp][p][q]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += -_flux_jacobi_grad_variable_ee[_qp][p][q]*_grad_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_ee[_qp][p][q]*_grad_test_neighbor[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::NeighborNeighbor:
		r = -_flux_jacobi_variable_en[_qp][p][q]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += -_flux_jacobi_grad_variable_en[_qp][p][q]*_grad_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_en[_qp][p][q]*_grad_test_neighbor[_i][_qp]*_phi_neighbor[_j][_qp];
		break;
	}

	return r;
}
