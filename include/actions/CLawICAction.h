
#pragma once

#include "Action.h"

class CLawICAction : public Action
{
public:
  CLawICAction(const std::string & name, InputParameters params);
  virtual void act();

protected:
};

template<>
InputParameters validParams<CLawICAction>();
