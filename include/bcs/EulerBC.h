
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
	  MaterialProperty<std::vector<RealVectorValue> > &_invis_term;
	  MaterialProperty<std::vector<RealVectorValue> > &_invis_term_neighbor;
	  MaterialProperty<Real > &_flux_diff;

	  virtual Real computeQpResidual();

};
