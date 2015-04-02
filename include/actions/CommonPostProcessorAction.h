
#pragma once

#include "Action.h"


class CommonPostProcessorAction: public Action
{
public:
  CommonPostProcessorAction(const std::string & name, InputParameters params);
  virtual void act();

private:
  InputParameters _action_params;
};

template<>
InputParameters validParams<CommonPostProcessorAction>();

