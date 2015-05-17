
#include "CLawBoundaryCondition.h"

template<>
InputParameters validParams<CLawBoundaryCondition>()
{
	InputParameters params = validParams<IntegratedBC>();
	params.addRequiredParam<int>("component", "BC component");
	return params;
}

CLawBoundaryCondition::CLawBoundaryCondition(const std::string & name, InputParameters parameters):
		IntegratedBC(name, parameters),
		_boundary(getMaterialProperty<CLawBoundaryMaterialData>("bnd_material_data")),
		_eq(getParam<int>("component"))
{
}

Real CLawBoundaryCondition::computeQpResidual()
{
	return _boundary[_qp]._flux[_eq] * _test[_i][_qp] +
			_boundary[_qp]._lift[_eq]*_grad_test[_i][_qp];
}

Real CLawBoundaryCondition::computeQpJacobian()
{
	int p(_eq), q(_eq);
	return computeQpJacobian(p, q);
}

Real CLawBoundaryCondition::computeQpOffDiagJacobian(unsigned int jvar)
{
	int p(_eq), q(jvar);
	return computeQpJacobian(p, q);
}

Real CLawBoundaryCondition::computeQpJacobian(int p, int q)
{
	Real r(0);

	r = _boundary[_qp]._flux_jacobi_variable[p][q]*_phi[_j][_qp]*_test[_i][_qp];
	r += _boundary[_qp]._flux_jacobi_grad_variable[p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
	r += _boundary[_qp]._lift_jacobi_variable[p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];

	return r;
}
