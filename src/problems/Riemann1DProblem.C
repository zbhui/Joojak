
#include "Riemann1DProblem.h"
#include "CLawBoundaryMaterial.h"
#include "boost/assign.hpp"
#include "boost/typeof/typeof.hpp"
#include <functional>
using namespace boost::assign;

template<>
InputParameters validParams<Riemann1DProblem>()
{
  InputParameters params = validParams<EulerProblem>();
  MooseEnum sub_types("sod lax shu blast");
  params.addRequiredParam<MooseEnum>("sub_type", sub_types, "1D激波管");
  return params;
}

Riemann1DProblem::Riemann1DProblem(const std::string & name, InputParameters params) :
	EulerProblem(name, params),
	_sub_type(getParam<MooseEnum>("sub_type"))
{
	if(_sub_type == "sod")
	{
		_initial_condition[0] += 1,0,0,0, 1;
		_initial_condition[1] += 0.125,0,0,0, 0.1;

		std::vector<Point> points;
		points += Point(0), Point(0.5), Point(1);
		_initial[MeshTools::BoundingBox(points[0], points[1])] = _initial_condition[0];
		_initial[MeshTools::BoundingBox(points[1], points[2])] = _initial_condition[1] ;
	}
	if(_sub_type == "lax")
	{
		_initial_condition[0] += 0.455,0.698,0,0, 3.528;
		_initial_condition[1] += 0.5,0,0,0, 0.571;

		std::vector<Point> points;
		points += Point(0), Point(0.5), Point(1);
		_initial[MeshTools::BoundingBox(points[0], points[1])] = _initial_condition[0];
		_initial[MeshTools::BoundingBox(points[1], points[2])] = _initial_condition[1] ;
	}
	if(_sub_type == "blast")
	{
		_initial_condition[0] += 1, 0, 0, 0, 1000;
		_initial_condition[1] += 1, 0, 0, 0, 0.1;
		_initial_condition[2] += 1, 0, 0, 0, 100;

		std::vector<Point> points;
		points += Point(0), Point(0.1), Point(0.9), Point(1);
		_initial[MeshTools::BoundingBox(points[0], points[1])] = _initial_condition[0];
		_initial[MeshTools::BoundingBox(points[1], points[2])] = _initial_condition[1] ;
		_initial[MeshTools::BoundingBox(points[2], points[3])] = _initial_condition[1] ;
	}
	if(_sub_type == "shu")
	{
		_initial_condition[0] += 3.857143,2.629369,0,0, 10.3333333;
		_initial_condition[1] += 100000,0,0,0, 1;

		std::vector<Point> points;
		points += Point(-5), Point(-4), Point(5);
		_initial[MeshTools::BoundingBox(points[0], points[1])] = _initial_condition[0];
		_initial[MeshTools::BoundingBox(points[1], points[2])] = _initial_condition[1] ;
	}
}

vector<Real> & Riemann1DProblem::valueAtPoint(const Point &p)
{
	for(BOOST_AUTO(it, _initial.begin()); it != _initial.end(); ++it)
	{
		if(it->first.contains_point(p))
			return it->second;
	}
	mooseError("初始条件越界了");
}

Real Riemann1DProblem::density(Real t, const Point &p)
{
	vector<Real> value= valueAtPoint(p);
	if(_sub_type == "shu")
	{
		if(p(0) > -4)
			value[0] = 1+0.2*sin(5*p(0));
	}
	return value[0];
}

Real Riemann1DProblem::momentumX(Real t, const Point &p)
{
	vector<Real> value = valueAtPoint(p);
	return value[0]*value[1];
}

Real Riemann1DProblem::momentumY(Real t, const Point &p)
{
	return 0.0;
}

Real Riemann1DProblem::momentumZ(Real t, const Point &p)
{
	return 0.0;
}

Real Riemann1DProblem::energyTotal(Real t, const Point &p)
{
	RealVectorValue momentum(momentumX(t, p), momentumY(t, p), momentumZ(t, p));
	Real rho = density(t, p);
	Real pre = pressure(t, p);

	return pre/(_gamma-1) +0.5*momentum.size_sq()/rho;
}

Real Riemann1DProblem::pressure(Real t, const Point &p)
{
	vector<Real> value = valueAtPoint(p);
	return value[4];
}
void Riemann1DProblem::computeBoundaryFlux(Real* flux, RealVectorValue* lift, Real* ul, RealGradient* dul, CLawBoundaryMaterial& bnd)
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
		ur[eq] = valueExact(time(), q_point, eq);

	if(_sub_type == "blast")
		EulerProblem::wall(ur, ul, normal);

	computeFaceFlux(flux, lift, ul, ur, dul, dur, normal, penalty);
}

Real Riemann1DProblem::initialCondition(const Point& point, int eq)
{
	return valueExact(0, point, eq);
}

Real Riemann1DProblem::valueExact(Real t, const Point& p, int eq)
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
