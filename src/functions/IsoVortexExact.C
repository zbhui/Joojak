/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "IsoVortexExact.h"

template<>
InputParameters validParams<IsoVortexExact>()
{
  InputParameters params = validParams<Function>();
  params += validParams<IsoVortexBase>();
  return params;
}

IsoVortexExact::IsoVortexExact(const std::string & name, InputParameters parameters) :
    Function(name, parameters),
    IsoVortexBase(name, parameters)
{
//	std::string var_name = _var.name();
//
//	if(var_name == "rho")
//		_eq = 0;
//	if(var_name == "momentum_x")
//		_eq = 1;
//	if(var_name == "momentum_y")
//		_eq = 2;
//	if(var_name == "momentum_z")
//		_eq = 3;
//	if(var_name == "rhoe")
//		_eq = 4;
}

Real IsoVortexExact::value(Real t, const Point & p)
{
	return IsoVortexBase::value(t, p, 0);
}
