
#pragma once

#include "SideUserObject.h"

class CFDForceUserObject : public SideUserObject
{
public:
	CFDForceUserObject(const std::string & name, InputParameters parameters);

	virtual void initialize(){};
	virtual void finalize();
	virtual void execute();
	virtual void threadJoin(const UserObject & uo){};
};


template<>
InputParameters validParams<CFDForceUserObject>();
