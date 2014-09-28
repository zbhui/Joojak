
#include "NearestWallDistance.h"

template<>
InputParameters validParams<NearestWallDistance>()
{
  InputParameters params = validParams<AuxKernel>();

  params.addRequiredCoupledVar("potential", "势函数");

  return params;
}

NearestWallDistance::NearestWallDistance(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
	_psi(coupledValue("potential")),
	_grad_psi(coupledGradient("potential"))
{}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per elemenet or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */
Real
NearestWallDistance::computeValue()
{
	return -_grad_psi[_qp].size()+sqrt(_grad_psi[_qp].size_sq()+2*_psi[_qp]);
}
