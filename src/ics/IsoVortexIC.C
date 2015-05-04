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

#include "IsoVortexIC.h"

template<>
InputParameters validParams<IsoVortexIC>()
{
  InputParameters params = validParams<CFDInitialCondition>();
  params += validParams<IsoVortexBase>();
  return params;
}

IsoVortexIC::IsoVortexIC(const std::string & name, InputParameters parameters) :
    CFDInitialCondition(name, parameters),
    IsoVortexBase(name, parameters)
{}

Real IsoVortexIC::density(const Point &p)
{
	return IsoVortexBase::density(0, p);
}

Real IsoVortexIC::momentumX(const Point &p)
{
	return IsoVortexBase::x_momentum(0, p);
}

Real IsoVortexIC::momentumY(const Point &p)
{
	return IsoVortexBase::y_momentum(0, p);
}

Real IsoVortexIC::momentumZ(const Point &p)
{
	return IsoVortexBase::z_momentum(0, p);
}

Real IsoVortexIC::energyTotal(const Point &p)
{
	return IsoVortexBase::total_energy(0, p);
}
