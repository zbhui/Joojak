
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawProblem>()
{
  InputParameters params = validParams<FEProblem>();
  params.addCoupledVar("aux_variables", "耦合变量");
  return params;
}

CLawProblem::CLawProblem(const std::string & name, InputParameters params) :
    FEProblem(name, params),
	_variables(_nl.getVariableNames()),
	_n_equations(_variables.size())
{
	std::cout << "配置文件：" << _app.getInputFileName() << std::endl;
	std::cout << params <<std::endl;
	std::cout << _n_equations <<std::endl;
//    _mesh.printInfo();
}



void CLawProblem::init()
{
	FEProblem::init();
	_variables = _nl.getVariableNames();
	_n_equations = (_variables.size());
	std::cout << _n_equations <<std::endl;
}
void CLawProblem::computeFaceFlux(Real* flux, RealVectorValue* lift, Real* ul, Real* ur, RealGradient* dul, RealGradient* dur, Point& normal, Real penalty)
{
	RealVectorValue ifl[10], ifr[10], vfl[10], vfr[10];
	computeLift(lift, ul, ur, normal);

	viscousTerm(vfl, ul, dul);
	viscousTerm(vfr, ur, dur);

	fluxRiemann(flux, ul, ur, normal);

	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux[eq] -= 0.5*((vfl[eq]+vfr[eq])-penalty*lift[eq])*normal;
	}

}

void CLawProblem::computeLift(RealVectorValue *lift, Real *ul, Real *ur, Point &normal)
{
	RealGradient duh[10];
	Real uh[10];
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		duh[eq] = (ul[eq]-ur[eq])/2.*normal;
		uh[eq] = (ul[eq]+ur[eq])/2;
	}
	viscousTerm(lift, uh, duh);
}

void CLawProblem::computeBoundaryFlux(Real* flux, RealVectorValue* lift, Real* ul, RealGradient* dul, Point& normal, Real penalty, std::string bc_type)
{
	Real ur[10];
	RealGradient dur[10];
	RealVectorValue ifl[10], ifr[10], vfl[10], vfr[10];

	for (int eq = 0; eq < _n_equations; ++eq)
		dur[eq] = dul[eq];

	boundaryCondition(ur, ul, normal, bc_type);

	computeFaceFlux(flux, lift, ul, ur, dul, dur, normal, penalty);
}

void CLawProblem::computeCellFlux(RealGradient* flux, Real* uh, RealGradient* duh)
{
	RealVectorValue inv_term[10], vis_term[10], source_term[10];
	inviscousTerm(inv_term, uh);
	viscousTerm(vis_term, uh, duh);
//	sourceTerm(source_term, uh, duh);
	for (int eq = 0; eq < _n_equations; ++eq)
		flux[eq] = inv_term[eq] - vis_term[eq];
}

void CLawProblem::inviscousTerm(RealVectorValue* inviscous_term, Real* uh)
{
	mooseError("CLawProblem::inviscousTerm，需要子类填充.");

}

void CLawProblem::sourceTerm(RealVectorValue* source_term, Real* uh, RealGradient* duh)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		source_term[eq](0) = 0;
		source_term[eq](1) = 0;
		source_term[eq](2) = 0;
	}
}

void CLawProblem::viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient* duh)
{
	mooseError("CLawProblem::viscousTerm，需要子类填充.");
}

void CLawProblem::fluxRiemann(Real* flux, Real* ul, Real* ur, Point& normal)
{
	mooseError("CLawProblem::fluxRiemann，需要子类填充.");
}

void CLawProblem::boundaryCondition(Real* ur, Real* ul, Point& normal, std::string bc_type)
{
	mooseError("CLawProblem::boundaryCondition，需要子类填充.");
}

Real CLawProblem::initialCondition(int eq)
{
	mooseError("CLawProblem::initialCondition，需要子类填充.");
}

