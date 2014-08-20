
#pragma once

#include "CFDBC.h"

class EulerBC;

template<>
InputParameters validParams<EulerBC>();

class EulerBC :
public CFDBC
{
public:
	  EulerBC(const std::string & name, InputParameters params);

protected:
	  MaterialProperty<std::vector<Real> > &_flux;
	  virtual Real computeQpResidual();

};
