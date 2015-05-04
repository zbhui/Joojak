#include "IsoVortexElementL2Error.h"
#include "EulerProblem.h"

template<>
InputParameters validParams<IsoVortexElementL2Error>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params += validParams<IsoVortexBase>();
  return params;
}

IsoVortexElementL2Error::IsoVortexElementL2Error(const std::string & name, InputParameters parameters) :
	ElementIntegralPostprocessor(name, parameters),
    IsoVortexBase(name, parameters)
{
	for (int eq = 0; eq < _euler_problem._n_equations; ++eq)
	{
		MooseVariable &val = _euler_problem.getVariable(_tid, "rho");
		_uh.push_back(_is_implicit ? &val.sln() : &val.slnOld());
	}
}

Real IsoVortexElementL2Error::getValue()
{
  return std::sqrt(ElementIntegralPostprocessor::getValue());
}

Real IsoVortexElementL2Error::computeQpIntegral()
{
	Real uh[5];
	for (int eq = 0; eq < _euler_problem._n_equations; ++eq)
		uh[eq] = (*_uh[eq])[_qp];

	Real err = uh[0] - value(_t, _q_point[_qp] , 0);
	return err*err;
}
