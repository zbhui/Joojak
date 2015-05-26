
#include "FluxJumpIndicator.h"
#include "CLawProblem.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<FluxJumpIndicator>()
{
  InputParameters params = validParams<InternalSideIndicator>();
  return params;
}


FluxJumpIndicator::FluxJumpIndicator(const std::string & name, InputParameters parameters) :
	InternalSideIndicator(name, parameters),
	_claw_problem(static_cast<CLawProblem&>(*parameters.get<FEProblem *>("_fe_problem"))),
	_cfd_problem(static_cast<CFDProblem&>(*parameters.get<FEProblem *>("_fe_problem"))),
	_nl(_claw_problem.getNonlinearSystem()),
	_tid(parameters.get<THREAD_ID>("_tid")),
	_variables(_nl.getVariableNames()),
	_var_order(_claw_problem.getVariable(_tid, _variables[0]).order()),
	_is_implicit(true)
{
	_aux_variables = _claw_problem._aux_variables;
	for(int i = 0; i < _aux_variables.size(); ++i)
		_variables.push_back(_aux_variables[i]);
	_n_variables = _variables.size();

	for (int ivar = 0; ivar < _n_variables; ++ivar)
	{
		MooseVariable &val = _claw_problem.getVariable(_tid, _variables[ivar]);
		_uh.push_back(_is_implicit ? &val.sln() : &val.slnOld());
		_uh_neighbor.push_back(_is_implicit ? &val.slnNeighbor() : &val.slnOldNeighbor());
		_grad_uh.push_back(_is_implicit ? &val.gradSln(): &val.gradSlnOld());
		_grad_uh_neighbor.push_back(_is_implicit ? &val.gradSlnNeighbor(): &val.gradSlnOldNeighbor());
	}
}

void FluxJumpIndicator::computeIndicator()
{
  Real sum = 0;

  for (_qp=0; _qp<_qrule->n_points(); _qp++)
    sum += _JxW[_qp]*_coord[_qp]*computeQpIntegral();

  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

    _solution.add(_field_var.nodalDofIndex(), sum/_current_side_elem->volume());
    _solution.add(_field_var.nodalDofIndexNeighbor(), sum/_current_side_elem->volume());

  }
}

Real FluxJumpIndicator::computeQpIntegral()
{
	Real ul[10], ur[10], uh[10];
	RealVectorValue ifl[10], ifr[10];
	for (int eq = 0; eq < _n_variables; ++eq)
	{
		ul[eq] =  (*_uh[eq])[_qp];
		ur[eq] =  (*_uh_neighbor[eq])[_qp];
		uh[eq] = (ul[eq] + ur[eq])/2.;
	}

	_claw_problem.inviscousTerm(ifl, ul);
	_claw_problem.inviscousTerm(ifr, ur);

	Real pressure = _cfd_problem.pressure(uh);
	Real pressure_new(0);
	Real ds = 1e-08;
	Real indicator(0);
	for (int eq = 0; eq < _claw_problem._n_equations; ++eq)
	{
		Real flux_jump= (ifl[eq] - ifr[eq]) * _normals[_qp];

		uh[eq] += ds;
		pressure_new = _cfd_problem.pressure(uh);

		indicator += (flux_jump/pressure * (pressure_new - pressure)/ds);
	}

	return (indicator);
}

void FluxJumpIndicator::finalize()
{
  unsigned int n_flux_faces = 0;

  if (_scale_by_flux_faces)
  {
    // Figure out the total number of sides contributing to the error.
    // We'll scale by this so boundary elements are less penalized
    for (unsigned int side=0; side<_current_elem->n_sides(); side++)
      if (_current_elem->neighbor(side) != NULL)
        n_flux_faces++;
  }
  else
    n_flux_faces = 1;

  // The 0 is because CONSTANT MONOMIALS only have one coefficient per element...
  Real value = _field_var.nodalSln()[0];

  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    _solution.set(_field_var.nodalDofIndex(), std::fabs(value)/(Real)n_flux_faces);
  }
}

