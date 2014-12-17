
#pragma once

#include "Action.h"
#include "libmesh/fe.h"

class CFDMaterialAction;

template<>
InputParameters validParams<CFDMaterialAction>();

class CFDMaterialAction : public Action
{
public:
	CFDMaterialAction(const std::string & name, InputParameters params);
  virtual void act();

protected:

};
