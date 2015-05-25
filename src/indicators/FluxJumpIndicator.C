
#include "FluxJumpIndicator.h"

template<>
InputParameters validParams<FluxJumpIndicator>()
{
  InputParameters params = validParams<InternalSideIndicator>();
  return params;
}


FluxJumpIndicator::FluxJumpIndicator(const std::string & name, InputParameters parameters) :
	InternalSideIndicator(name, parameters)
{
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
  Real jump = fabs(_u[_qp] - _u_neighbor[_qp]);

  return jump;
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

