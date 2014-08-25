
#pragma once

#include "IntegratedBC.h"
#include "MooseEnum.h"

class CFDBC;

template<>
InputParameters validParams<CFDBC>();

class CFDBC :
public IntegratedBC
{
public:
	CFDBC(const std::string & name, InputParameters params);

protected:

	 MooseEnum _bc_type;
};


