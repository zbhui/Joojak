
#include "RatioTimeStepper.h"

//Moose includes

template<>
InputParameters validParams<RatioTimeStepper>()
{
  InputParameters params = validParams<TimeStepper>();
  params.addParam<Real>("dt", 1., "The initial time step size.");
  params.addParam<Real>("ratio", 1, "The ratio used to calculate the next timestep");
  params.addParam<Real>("max_dt", 1e09, "The largest timestep we will allow");
  params.addParam<int> ("step", 1, "时间步长更新频度");

  return params;
}

RatioTimeStepper::RatioTimeStepper(const std::string & name, InputParameters parameters) :
    TimeStepper(name, parameters),
    _ratio(getParam<Real>("ratio")),
    _max_dt(getParam<Real>("max_dt")),
    _step(getParam<int>("step"))
{
	mooseAssert(_step > 0, "更新频度为负数");
}

Real
RatioTimeStepper::computeInitialDT()
{
  return getParam<Real>("dt");
}

Real
RatioTimeStepper::computeDT()
{
  if (_ratio > 1.0)
	 if(_t_step % _step == 0)
		 return std::min(getCurrentDT() * _ratio, _max_dt);
	 else
		 return getCurrentDT();
  else
    return getCurrentDT();
}
