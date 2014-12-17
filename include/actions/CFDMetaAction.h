
#pragma once
#include "Action.h"

class CFDMetaAction : public Action
{
public:
	CFDMetaAction(const std::string& name, InputParameters params);

  virtual void act();

protected:
  const std::vector<BoundaryName> _boundary;
};

template<>
InputParameters validParams<CFDMetaAction>();



