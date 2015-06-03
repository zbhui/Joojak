
#pragma once

#include "EulerProblem.h"
#include "libmesh/mesh_tools.h"
using std::vector;

class Riemann1DProblem : public EulerProblem
{
public:
	virtual Real initialCondition(const Point & point, int eq);
	Riemann1DProblem(const std::string & name, InputParameters params);

	virtual void computeBoundaryFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul, CLawBoundaryMaterial &bnd);

	Real valueExact(Real t, const Point &p, int eq);
private:
	Real density(Real t, const Point &p);
	Real momentumX(Real t, const Point &p);
	Real momentumY(Real t, const Point &p);
	Real momentumZ(Real t, const Point &p);
	Real energyTotal(Real t, const Point &p);
	Real pressure(Real t, const Point &p);

	vector<Real> & valueAtPoint(const Point &p);
	vector<Real> _initial_condition[4];
	std::map<MeshTools::BoundingBox, vector<Real> > _initial;
	MooseEnum _sub_type;
};

template<>
InputParameters validParams<Riemann1DProblem>();
