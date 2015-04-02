
#pragma once

#include "Action.h"


class CommonPostProcessorAction: public Action
{
public:
  CommonPostProcessorAction(const std::string & name, InputParameters params);
  virtual void act();

private:
  void create(std::string type);
};

template<>
InputParameters validParams<CommonPostProcessorAction>();

