
#pragma once

#include "CFDBC.h"

class IsoVortexBC;

template<>
InputParameters validParams<IsoVortexBC>();

class IsoVortexBC :
public EulerBC
{
public:
	IsoVortexBC(const std::string & name, InputParameters params);
	virtual ~IsoVortexBC(){}


protected:
	virtual void valueAtRightFace(Real *ur);

private:

  Real density(Real t, const Point &p);
  Real x_momentum(Real t, const Point &p);
  Real y_momentum(Real t, const Point &p);
  Real z_momentum(Real t, const Point &p);
  Real total_energy(Real t, const Point &p);
};
