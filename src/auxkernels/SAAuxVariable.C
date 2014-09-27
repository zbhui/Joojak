
#include "SAAuxVariable.h"

template<>
InputParameters validParams<SAAuxVariable>()
{
  InputParameters params = validParams<AuxKernel>();
  params += validParams<SABase>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

SAAuxVariable::SAAuxVariable(const std::string & name, InputParameters parameters) :
    AuxKernel(name, parameters),
    SABase(name, parameters)
{
}

Real SAAuxVariable::computeValue()
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
		return eddyViscosity(uh);

	mooseError(var_name << "辅助变量名不存在");
	return 0.;

}

void SAAuxVariable::valueAtCellPoint(Real *uh)
{
	size_t n_equation = coupledComponents("variables");
	for (size_t eq = 0; eq < n_equation; ++eq)
	{
		uh[eq] = coupledValue("variables", eq)[_qp];
	}
}
