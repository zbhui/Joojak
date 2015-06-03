
#include "Riemann2DProblem.h"
#include "CLawBoundaryMaterial.h"
#include "boost/assign.hpp"
using namespace boost::assign;

template<>
InputParameters validParams<Riemann2DProblem>()
{
  InputParameters params = validParams<EulerProblem>();

  return params;
}

Riemann2DProblem::Riemann2DProblem(const std::string & name, InputParameters params) :
	EulerProblem(name, params)
{
//	_initial_condition[0] = {0.5313,0,0,0,0.4};
//	_initial_condition[1] = {1,0.7276,0,0,1};
//	_initial_condition[2] = {0.8,0,0,0,1};
//	_initial_condition[3] = {1,0,0.7276,0,1};
	_initial_condition[0] += 0.5313,0,0,0,0.4;
	_initial_condition[1] += 1,0.7276,0,0,1;
	_initial_condition[2] += 0.8,0,0,0,1;
	_initial_condition[3] += 1,0,0.7276,0,1;
}

Real Riemann2DProblem::density(Real t, const Point &p)
{
	return _initial_condition[pointLocator(p)][0];
}

Real Riemann2DProblem::momentumX(Real t, const Point &p)
{
	Real u = _initial_condition[pointLocator(p)][1];
	Real rho = density(t, p);

	return rho*u;
}

Real Riemann2DProblem::momentumY(Real t, const Point &p)
{
	Real u = _initial_condition[pointLocator(p)][2];
	Real rho = density(t, p);

	return rho*u;
}

Real Riemann2DProblem::momentumZ(Real t, const Point &p)
{
	Real u = _initial_condition[pointLocator(p)][3];
	Real rho = density(t, p);

	return rho*u;
}

Real Riemann2DProblem::energyTotal(Real t, const Point &p)
{
	RealVectorValue momentum(momentumX(t, p), momentumY(t, p), momentumZ(t, p));
	Real rho = density(t, p);
	Real pre = pressure(t, p);

	return pre/(_gamma-1) +0.5*momentum.size_sq()/rho;
}

Real Riemann2DProblem::pressure(Real t, const Point &p)
{
	return _initial_condition[pointLocator(p)][4];
}
void Riemann2DProblem::computeBoundaryFlux(Real* flux, RealVectorValue* lift, Real* ul, RealGradient* dul, CLawBoundaryMaterial& bnd)
{
	Point normal = bnd.normals()[_qp];
	Real penalty = bnd.penalty();
	RealGradient dur[10];
	Real ur[10];
	for (int eq = 0; eq < _n_equations; ++eq)
		dur[eq] = dul[eq];

	std::string bc_type = bnd.getBCType();
	Point q_point = bnd.qpoints()[_qp];
	for (int eq = 0; eq < _n_equations; ++eq)
//		ur[eq] = valueExact(time(), q_point, eq);
		ur[eq] = ul[eq];

	computeFaceFlux(flux, lift, ul, ur, dul, dur, normal, penalty);
}

Real Riemann2DProblem::initialCondition(const Point& point, int eq)
{
	return valueExact(0, point, eq);
}

Real Riemann2DProblem::valueExact(Real t, const Point& p, int eq)
{
	switch (eq) {
	case 0:
		return density(t, p);
		break;
	case 1:
		return momentumX(t, p);
		break;
	case 2:
		return momentumY(t, p);
		break;
	case 3:
		return momentumZ(t, p);
	case 4:
		return energyTotal(t, p);
		break;
	default:
		return 0.0;
		mooseError("不可用的分量" << eq);
		break;
	}
}

int Riemann2DProblem::pointLocator(const Point& p)
{
	Real x = p(0);
	Real y = p(1);
    if (x > 0.5 && y > 0.5)
		return 0;
    else if (x < 0.5 && y > 0.5)
    	return 1;
    else if (x < 0.5 && y < 0.5)
    	return 2;
    else
    	return 3;
}
