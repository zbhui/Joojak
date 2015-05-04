
#pragma once

// MOOSE Includes
#include "CFDInitialCondition.h"
#include "IsoVortexBase.h"

// Forward Declarations
class IsoVortexIC;

template<>
InputParameters validParams<IsoVortexIC>();

/**
 * 等熵涡初始条件
 */
class IsoVortexIC :
public CFDInitialCondition,
public IsoVortexBase
{
public:
  IsoVortexIC(const std::string & name, InputParameters parameters);
protected:
  virtual Real density(const Point &p);
  virtual Real momentumX(const Point &p);
  virtual Real momentumY(const Point &p);
  virtual Real momentumZ(const Point &p);
  virtual Real energyTotal(const Point &p);
};
