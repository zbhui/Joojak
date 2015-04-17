

#pragma once

#include "InputParameters.h"
#include "Assembly.h"
using std::vector;

class CLawProblem;
class NonlinearSystem;

class CLawCellMaterialData
{
public:
	void update(CLawProblem & claw_problem);
	void setProblem(CLawProblem & claw_problem, Real ds);

private:
	void computeQpValue(RealVectorValue *flux_term);
	CLawProblem * _claw_problem;
	int _n_equations;
	Real _ds;
public:
	Real uh[10];
	RealGradient duh[10];
	RealVectorValue _flux_term[10];
	RealVectorValue _flux_jacobi_variable[10][10];
	RealTensorValue _flux_jacobi_grad_variable[10][10];
};
