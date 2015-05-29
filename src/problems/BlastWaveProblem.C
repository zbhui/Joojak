
#include "BlastWaveProblem.h"
#include "CLawBoundaryMaterial.h"

template<>
InputParameters validParams<BlastWaveProblem>()
{
  InputParameters params = validParams<EulerProblem>();

  return params;
}

BlastWaveProblem::BlastWaveProblem(const std::string & name, InputParameters params) :
	EulerProblem(name, params)
{
}

Real BlastWaveProblem::density(Real t, const Point &p)
{
	return 1;
}

Real BlastWaveProblem::momentumX(Real t, const Point &p)
{
	return 0;
}

Real BlastWaveProblem::momentumY(Real t, const Point &p)
{
	return 0.0;
}

Real BlastWaveProblem::momentumZ(Real t, const Point &p)
{
	return 0.0;
}

Real BlastWaveProblem::energyTotal(Real t, const Point &p)
{
	RealVectorValue momentum(momentumX(t, p), momentumY(t, p), momentumZ(t, p));
	Real rho = density(t, p);
	Real pre = pressure(t, p);

	return pre/(_gamma-1) +0.5*momentum.size_sq()/rho;
}

Real BlastWaveProblem::pressure(Real t, const Point &p)
{
	Real x = p(0);
    if (x < 10)
		return 1000;
    else if(x > 90)
    	return 100;
    else
    	return 0.01;
}

Real BlastWaveProblem::initialCondition(const Point& point, int eq)
{
	return valueExact(0, point, eq);
}

Real BlastWaveProblem::valueExact(Real t, const Point& p, int eq)
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
