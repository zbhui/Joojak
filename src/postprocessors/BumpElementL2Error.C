#include "BumpElementL2Error.h"

template<>
InputParameters validParams<BumpElementL2Error>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params += validParams<EulerBase>();
  params.addRequiredParam<PostprocessorName>("area", "归一化面积");
  params.addRequiredCoupledVar("variables", "守恒变量");
  return params;
}

BumpElementL2Error::BumpElementL2Error(const std::string & name, InputParameters parameters) :
	ElementIntegralPostprocessor(name, parameters),
    EulerBase(name, parameters),
    _area(getPostprocessorValue("area"))
{
}

Real BumpElementL2Error::getValue()
{
  return std::sqrt(ElementIntegralPostprocessor::getValue()/_area);
}

Real BumpElementL2Error::computeQpIntegral()
{
	Real uh[5];
	int _n_equations = coupledComponents("variables");
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		uh[eq] = coupledValue("variables", eq)[_qp];
	}

	Real pre = pressure(uh);
	Real entropy = pre/std::pow(uh[0] ,_gamma);
	Real pre_inf = 1/_gamma/_mach/_mach;
	Real rho_inf = 1;
	Real entropy_inf = pre_inf/std::pow(rho_inf, _gamma);
	Real err = 1-entropy/entropy_inf;
	return err*err;
}
