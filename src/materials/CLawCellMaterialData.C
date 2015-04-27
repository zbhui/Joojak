
#include "CLawCellMaterialData.h"
#include "CLawProblem.h"

void CLawCellMaterialData::update(CLawProblem & claw_problem)
{
	RealVectorValue flux_term_new[10], vis_term[10];

	computeQpValue(_flux_term);

	for (int q = 0; q < _n_equations; ++q)
	{
		uh[q] += _ds;
		computeQpValue(flux_term_new);
		for (int p = 0; p < _n_equations; ++p)
			_flux_jacobi_variable[p][q] = (flux_term_new[p] - _flux_term[p])/_ds;

		uh[q] -= _ds;

		for (int beta = 0; beta < 3; ++beta)
		for (int q = 0; q < _n_equations; ++q)
		{
			duh[q](beta) += _ds;

			computeQpValue(flux_term_new);
			for (int alpha = 0; alpha< 3; ++alpha)
			for (int p = 0; p < _n_equations; ++p)
			{
				_flux_jacobi_grad_variable[p][q](alpha, beta) = (flux_term_new[p](alpha) - _flux_term[p](alpha))/_ds;
			}

			duh[q](beta) -= _ds;
		}
	}
}

void CLawCellMaterialData::computeQpValue(RealVectorValue *flux_term)
{
//	RealVectorValue inv_term[10], vis_term[10];
//	_claw_problem->inviscousTerm(inv_term, uh);
//	_claw_problem->viscousTerm(vis_term, uh, duh);
//	for (int eq = 0; eq < _n_equations; ++eq)
//		flux_term[eq] = inv_term[eq] - vis_term[eq];

	_claw_problem->computeCellFlux(flux_term, uh, duh);
}

void CLawCellMaterialData::setProblem(CLawProblem& claw_problem, Real ds)
{
	_claw_problem = & claw_problem;
	_n_equations = claw_problem._n_equations;
	_ds = ds;
}
