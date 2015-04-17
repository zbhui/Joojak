
#include "CLawBoundaryMaterial.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawBoundaryMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredCoupledVar("variables", "守恒变量");

//  MooseEnum bc_types;//("isothermal_wall adiabatic_wall far_field symmetric pressure_out none", "none");  // 边界条件的类型，可以增加
  params.addRequiredParam<std::string>("bc_type", "边界条件");

  return params;
}

CLawBoundaryMaterial::CLawBoundaryMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		CLawInterface(parameters),
		_bc_type(getParam<std::string>("bc_type")),
		_flux(declareProperty<std::vector<Real> >("flux")),
		_flux_jacobi_variable(declareProperty<std::vector<std::vector<Real> > >("flux_jacobi_variable")),
		_flux_jacobi_grad_variable(declareProperty<std::vector<std::vector<RealGradient> > >("flux_jacobi_grad_variable")),

		_lift(declareProperty<std::vector<RealVectorValue> >("lift")),
		_lift_jacobi_variable(declareProperty<std::vector<std::vector<RealVectorValue> > >("lift_jacobi_variable_ee"))
{
	_n_equations = coupledComponents("variables");
	for (size_t eq = 0; eq < _n_equations; ++eq)
	{
		MooseVariable &val = *getVar("variables", eq);
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
	computeQpRightValue(ur, dur, ul, dul);

	liftOperator(&_lift[_qp][0], ul, ur);
	fluxTerm(&_flux[_qp][0], ul, ur, dul, dur);

	Real _ds = 1E-08;
	for (int q = 0; q < _n_equations; ++q)
	{
		ul[q] += _ds;
		computeQpRightValue(ur_new, dur_new, ul, dul);
		liftOperator(lift_new, ul, ur_new);
		fluxTerm(flux_new, ul, ur_new, dul, dur_new);
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
		computeQpRightValue(ur_new, dur_new, ul, dul);
		fluxTerm(flux_new, ul, ur_new, dul, dur_new);
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

void CLawBoundaryMaterial::computeQpRightValue(Real *ur, RealGradient *dur, Real *ul, RealGradient *dul)
{
	for (int q = 0; q < _n_equations; ++q)
		dur[q] = dul[q];

	Point normal = _normals[_qp];
	_claw_problem.boundaryCondition(ur, ul, normal, _bc_type);

}

void CLawBoundaryMaterial::fluxTerm(Real *flux, Real* ul, Real* ur, RealGradient *dul, RealGradient *dur)
{
	RealVectorValue ifl[5], ifr[5], vfl[5], vfr[5];

	convertionTerm(ifl, ul);
	convertionTerm(ifr, ur);
	diffusionTerm(vfl, ul, dul);
	diffusionTerm(vfr, ur, dur);

	Real lam = 1;//(maxEigenValue(ul, _normals[_qp]) + maxEigenValue(ur, _normals[_qp]))/2.;
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux[eq] = 0.5*(ifl[eq] + ifr[eq] - (vfl[eq]+vfl[eq]))*_normals[_qp] + lam*(ul[eq] - ur[eq]);
	}
}

void CLawBoundaryMaterial::liftOperator(RealVectorValue* lift, Real* ul, Real* ur)
{
	RealGradient duh[10];
	Real uh[10];
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		duh[eq] = (ul[eq]-ur[eq])/2.*_normals[_qp];
		uh[eq] = (ul[eq]+ur[eq])/2;
	}

	diffusionTerm(lift, uh, duh);
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
