
#include "NSAuxVariable.h"

template<>
InputParameters validParams<NSAuxVariable>()
{
  InputParameters params = validParams<AuxKernel>();
  params += validParams<NSBase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

NSAuxVariable::NSAuxVariable(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    NSBase(name, parameters)
{
	size_t n_equation = coupledComponents("variables");
	for (size_t i = 0; i < n_equation; ++i)
	{
		_uh.push_back(&coupledValue("variables", i));
	}

}

/**
 * Auxiliary Kernels override computeValue() instead of computeQpResidual().  Aux Variables
 * are calculated either one per elemenet or one per node depending on whether we declare
 * them as "Elemental (Constant Monomial)" or "Nodal (First Lagrange)".  No changes to the
 * source are necessary to switch from one type or the other.
 */
Real NSAuxVariable::computeValue()
{
	Real uh[10];
	valueAtCellPoint(uh);
	std::string var_name = _var.name();

	if(var_name == "pressure")
		return pressure(uh);
	if(var_name == "mach")
		return mach_local(uh);
	if(var_name == "velocity_x")
		return uh[1]/uh[0];
	if(var_name == "velocity_y")
		return uh[2]/uh[0];
	if(var_name == "velocity_z")
		return uh[3]/uh[0];
	if(var_name == "eddy_viscosity")
		return _reynolds * uh[5]/exp(uh[6]/uh[0]);

	mooseError(var_name << "辅助变量名不存在");
	return 0.;

}

void NSAuxVariable::valueAtCellPoint(Real *uh)
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		uh[i] = (*_uh[i])[_qp];
	}
}
