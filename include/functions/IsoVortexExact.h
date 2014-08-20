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

#ifndef IsoVortexExact_H
#define IsoVortexExact_H

#include "Function.h"
#include "IsoVortexBase.h"

class IsoVortexExact;

template<>
InputParameters validParams<IsoVortexExact>();

class IsoVortexExact :
public Function,
public IsoVortexBase
{
public:
  IsoVortexExact(const std::string & name, InputParameters parameters);

  Real value(Real t, const Point & p);

};

#endif //IsoVortexExact_H
