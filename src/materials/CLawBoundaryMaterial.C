
#include "CLawBoundaryMaterial.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawBoundaryMaterial>()
{
  InputParameters params = validParams<Material>();

  params.addRequiredParam<std::string>("bc_type", "边界条件");
  params.addParam<Real>("ds", 1.490116119384766e-08, "微扰量");
  params.addParam<Real>("sigma", 6, "通量罚值，默认值为6");
  params.addParam<Real>("epsilon", 1, "对称项罚值，可以取1, 0 , -1，分别对应SIP, IIP, NIP");

  return params;
}

CLawBoundaryMaterial::CLawBoundaryMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		CLawInterface(parameters),
		_bc_type(getParam<std::string>("bc_type")),
		_current_elem_volume(_assembly.elemVolume()),
		_neighbor_elem_volume(_assembly.neighborVolume()),
		_current_side_volume(_assembly.sideElemVolume()),
		_ds(getParam<Real>("ds")),
		_sigma(getParam<Real>("sigma")),
		_epsilon(getParam<Real>("epsilon")),
		_flux(declareProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable(declareProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable")),
		_flux_jacobi_grad_variable(declareProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable")),

		_lift(declareProperty<std::vector<RealVectorValue> >("lift")),
		_lift_jacobi_variable(declareProperty<std::vector<std::vector<RealVectorValue> > >("lift_jacobi_variable"))
{
	for (size_t eq = 0; eq < _n_equations; ++eq)
	{
		MooseVariable &val = getVariable(eq);
		_ul.push_back(_is_implicit ? &val.sln() : &val.slnOld());
		_grad_ul.push_back(_is_implicit ? &val.gradSln(): &val.gradSlnOld());
	}
}

void CLawBoundaryMaterial::computeQpProperties()
{
	if(!_bnd)
		mooseError("边界Material不在边界");

	resizeQpProperty();

	Real ul[10], ur[10], ur_new[10];
	RealGradient dul[10], dur[10], dur_new[10];
	Real flux_new[10];
	RealVectorValue lift_new[10];

	computeQpLeftValue(ul, dul);
	computeQpFlux(&_flux[_qp][0], &_lift[_qp][0], ul, dul);

	for (int q = 0; q < _n_equations; ++q)
	{
		ul[q] += _ds;
		computeQpFlux(flux_new, lift_new, ul, dul);
		for (int p = 0; p < _n_equations; ++p)
		{
			_flux_jacobi_variable[_qp][p][q] = (flux_new[p] - _flux[_qp][p])/_ds;
			_lift_jacobi_variable[_qp][p][q] = (lift_new[p] - _lift[_qp][p])/_ds;
		}
		ul[q] -= _ds;
	}

	for (int beta = 0; beta < 3; ++beta)
	for (int q = 0; q < _n_equations; ++q)
	{
		dul[q](beta) += _ds;
		computeQpFlux(flux_new, lift_new, ul, dul);
		for (int p = 0; p < _n_equations; ++p)
		{
			_flux_jacobi_grad_variable[_qp][p][q](beta) = (flux_new[p] - _flux[_qp][p])/_ds;
		}
		dul[q](beta) -= _ds;
	}

}

void CLawBoundaryMaterial::computeQpLeftValue(Real* ul, RealGradient *dul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		ul[eq] = (*_ul[eq])[_qp];
		dul[eq] = (*_grad_ul[eq])[_qp];
	}
}

void CLawBoundaryMaterial::computeQpFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul)
{
	Real h_face = (_current_elem_volume+_current_elem_volume)/_current_side_volume /2.;
	Real penalty = _sigma*_var_order*_var_order/h_face;
	Point normal = _normals[_qp];
	_claw_problem.computeBoundaryFlux(flux, lift, ul, dul, normal, penalty, _bc_type);
}

void CLawBoundaryMaterial::resizeQpProperty()
{
	_flux[_qp].resize(_n_equations);
	_flux_jacobi_variable[_qp].resize(_n_equations);
	_flux_jacobi_grad_variable[_qp].resize(_n_equations);
	_lift[_qp].resize(_n_equations);
	_lift_jacobi_variable[_qp].resize(_n_equations);
	for (int p = 0; p < _n_equations; ++p)
	{
		_flux_jacobi_variable[_qp][p].resize(_n_equations);
		_flux_jacobi_grad_variable[_qp][p].resize(_n_equations);
		_lift_jacobi_variable[_qp][p].resize(_n_equations);
	}
}
