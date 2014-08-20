
#pragma once

#include "IntegratedBC.h"
#include "CFDBase.h"
#include "MooseEnum.h"

class CFDBC;

template<>
InputParameters validParams<CFDBC>();

class CFDBC :
public IntegratedBC,
public CFDBase
{
public:
	CFDBC(const std::string & name, InputParameters params);

protected:

	 MooseEnum _bc_type;
	 int _eq;
};


