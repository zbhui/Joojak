
#pragma once

#include "TimeDerivative.h"

class EmptyTimeDerivative : public TimeDerivative
{
public:
	EmptyTimeDerivative(const std::string & name, InputParameters parameters);

  virtual void computeJacobian();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

};

template<>
InputParameters validParams<EmptyTimeDerivative>();
