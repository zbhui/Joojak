
#include "TestJumpIndicator.h"
#include "CLawProblem.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<TestJumpIndicator>()
{
  InputParameters params = validParams<InternalSideIndicator>();
  return params;
}


TestJumpIndicator::TestJumpIndicator(const std::string & name, InputParameters parameters) :
	InternalSideIndicator(name, parameters),
	_claw_problem(static_cast<CLawProblem&>(*parameters.get<FEProblem *>("_fe_problem"))),
	_cfd_problem(static_cast<CFDProblem&>(*parameters.get<FEProblem *>("_fe_problem"))),
	_nl(_claw_problem.getNonlinearSystem()),
	_tid(parameters.get<THREAD_ID>("_tid")),
	_variables(_nl.getVariableNames()),
	_n_variables(_variables.size()),
	_var_order(_claw_problem.getVariable(_tid, _variables[0]).order()),
	_is_implicit(true)
{
	for (int ivar = 0; ivar < _n_variables; ++ivar)
	{
		MooseVariable &val = _claw_problem.getVariable(_tid, _variables[ivar]);
		_uh.push_back(_is_implicit ? &val.sln() : &val.slnOld());
		_uh_neighbor.push_back(_is_implicit ? &val.slnNeighbor() : &val.slnOldNeighbor());
		_grad_uh.push_back(_is_implicit ? &val.gradSln(): &val.gradSlnOld());
		_grad_uh_neighbor.push_back(_is_implicit ? &val.gradSlnNeighbor(): &val.gradSlnOldNeighbor());
	}
}

void TestJumpIndicator::computeIndicator()
{
  Real sum = 0;

  for (_qp=0; _qp<_qrule->n_points(); _qp++)
    sum += _JxW[_qp]*_coord[_qp]*computeQpIntegral();

  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

    _solution.add(_field_var.nodalDofIndex(), fabs(sum)/_current_side_elem->volume());
    _solution.add(_field_var.nodalDofIndexNeighbor(), fabs(sum)/_current_side_elem->volume());

  }
}

Real TestJumpIndicator::computeQpIntegral()
{
    Real jump = (_u[_qp] - _u_neighbor[_qp])/(_u[_qp] + _u_neighbor[_qp])*2.0;
	return jump;
}

void TestJumpIndicator::finalize()
{

}

