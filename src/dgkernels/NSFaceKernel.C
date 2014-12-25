
#include "NSFaceKernel.h"

template<>
InputParameters validParams<NSFaceKernel>()
{
	  InputParameters params = validParams<DGKernel>();
	  params += validParams<NSBase>();
	  return params;
}
NSFaceKernel::NSFaceKernel(const std::string & name, InputParameters parameters):
		DGKernel(name, parameters),
		NSBase(name, parameters),
		_flux(getNeighborMaterialProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable_ee(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_ee")),
		_flux_jacobi_variable_en(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_en")),
		_flux_jacobi_variable_ne(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_ne")),
		_flux_jacobi_variable_nn(getNeighborMaterialProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable_nn")),
		_flux_jacobi_grad_variable_ee(getNeighborMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_ee")),
		_flux_jacobi_grad_variable_en(getNeighborMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_en")),
		_flux_jacobi_grad_variable_ne(getNeighborMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_ne")),
		_flux_jacobi_grad_variable_nn(getNeighborMaterialProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable_nn")),

		_penalty(getNeighborMaterialProperty<std::vector<RealVectorValue> >("penalty")),
		_penalty_neighbor(getNeighborMaterialProperty<std::vector<RealVectorValue> >("penalty_neighbor")),
		_penalty_jacobi_variable_ee(getNeighborMaterialProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_ee")),
		_penalty_jacobi_variable_en(getNeighborMaterialProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_en")),
		_penalty_jacobi_variable_ne(getNeighborMaterialProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_ne")),
		_penalty_jacobi_variable_nn(getNeighborMaterialProperty<std::vector<std::vector<RealVectorValue> > >("penalty_jacobi_variable_nn"))

{
	_eq = equationIndex(_var.name());
}

Real NSFaceKernel::computeQpResidual(Moose::DGResidualType type)
{
	Real CIP = computeCIP();
	Real flux = _flux[_qp][_eq];
	flux += CIP*(_penalty[_qp][_eq] + _penalty_neighbor[_qp][_eq])*_normals[_qp];
	if(type == Moose::Element)
	{
		return  flux * _test[_i][_qp] + _epsilon * _penalty[_qp][_eq]*_grad_test[_i][_qp];
	}
	if(type == Moose::Neighbor)
	{
		return -flux * _test_neighbor[_i][_qp] + _epsilon*_penalty_neighbor[_qp][_eq]*_grad_test_neighbor[_i][_qp];
	}
	mooseError("face flux error.");
	return 0.;
}

Real NSFaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
	Real r = 0;
	Real CIP = computeCIP();
	int p(_eq), q(_eq);
	switch (type)
	{
	case Moose::ElementElement:
		r = _flux_jacobi_variable_ee[_qp][p][q]*_phi[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_ee[_qp][p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
		r += CIP*(_penalty_jacobi_variable_ee[_qp][p][q] + _penalty_jacobi_variable_ne[_qp][p][q])*_normals[_qp]*_phi[_j][_qp]*_test[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_ee[_qp][p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::ElementNeighbor:
		r = _flux_jacobi_variable_en[_qp][p][q]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_en[_qp][p][q]*_grad_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += CIP*(_penalty_jacobi_variable_en[_qp][p][q] + _penalty_jacobi_variable_nn[_qp][p][q])*_normals[_qp]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_en[_qp][p][q]*_grad_test[_i][_qp]*_phi_neighbor[_j][_qp];
		break;

	case Moose::NeighborElement:
		r = _flux_jacobi_variable_ne[_qp][p][q]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += _flux_jacobi_grad_variable_ne[_qp][p][q]*_grad_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r -= CIP*(_penalty_jacobi_variable_ee[_qp][p][q] + _penalty_jacobi_variable_ne[_qp][p][q])*_normals[_qp]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_ne[_qp][p][q]*_grad_test_neighbor[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::NeighborNeighbor:
		r = _flux_jacobi_variable_nn[_qp][p][q]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += _flux_jacobi_grad_variable_nn[_qp][p][q]*_grad_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r -= CIP*(_penalty_jacobi_variable_en[_qp][p][q] + _penalty_jacobi_variable_nn[_qp][p][q])*_normals[_qp]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_nn[_qp][p][q]*_grad_test_neighbor[_i][_qp]*_phi_neighbor[_j][_qp];
		break;
	}

	return r;
}

Real NSFaceKernel::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
	Real r = 0;
	Real CIP = computeCIP();
	int p(_eq), q(jvar);
	switch (type)
	{
	case Moose::ElementElement:
		r = _flux_jacobi_variable_ee[_qp][p][q]*_phi[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_ee[_qp][p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
		r += CIP*(_penalty_jacobi_variable_ee[_qp][p][q] + _penalty_jacobi_variable_ne[_qp][p][q])*_normals[_qp]*_phi[_j][_qp]*_test[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_ee[_qp][p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::ElementNeighbor:
		r = _flux_jacobi_variable_en[_qp][p][q]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_en[_qp][p][q]*_grad_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += CIP*(_penalty_jacobi_variable_en[_qp][p][q] + _penalty_jacobi_variable_nn[_qp][p][q])*_normals[_qp]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_en[_qp][p][q]*_grad_test[_i][_qp]*_phi_neighbor[_j][_qp];
		break;

	case Moose::NeighborElement:
		r = _flux_jacobi_variable_ne[_qp][p][q]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += _flux_jacobi_grad_variable_ne[_qp][p][q]*_grad_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r -= CIP*(_penalty_jacobi_variable_ee[_qp][p][q] + _penalty_jacobi_variable_ne[_qp][p][q])*_normals[_qp]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_ne[_qp][p][q]*_grad_test_neighbor[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::NeighborNeighbor:
		r = _flux_jacobi_variable_nn[_qp][p][q]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += _flux_jacobi_grad_variable_nn[_qp][p][q]*_grad_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r -= CIP*(_penalty_jacobi_variable_en[_qp][p][q] + _penalty_jacobi_variable_nn[_qp][p][q])*_normals[_qp]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += _epsilon*_penalty_jacobi_variable_nn[_qp][p][q]*_grad_test_neighbor[_i][_qp]*_phi_neighbor[_j][_qp];
		break;
	}

	return r;
}

Real NSFaceKernel::computeCIP()
{
	const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
	const double h_elem = (_current_elem_volume+_neighbor_elem_volume)/_current_side_volume * 1./std::pow(elem_b_order, 2.)/2.;
	return _sigma/h_elem;
}
