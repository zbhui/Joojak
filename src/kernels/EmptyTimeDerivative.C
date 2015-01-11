
#include "EmptyTimeDerivative.h"

template<>
InputParameters validParams<EmptyTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}

EmptyTimeDerivative::EmptyTimeDerivative(const std::string & name, InputParameters parameters) :
	TimeDerivative(name, parameters)
{
}

Real
EmptyTimeDerivative::computeQpResidual()
{
  return 0;
}

Real
EmptyTimeDerivative::computeQpJacobian()
{
  return 0;
}

void
EmptyTimeDerivative::computeJacobian()
{
  if (_lumping)
  {
    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());

    for (_i = 0; _i < _test.size(); _i++)
      for (_j = 0; _j < _phi.size(); _j++)
        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
        {
          ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
        }
  }
  else
    TimeKernel::computeJacobian();
}
