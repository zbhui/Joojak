
#pragma once

#include "InternalSideIndicator.h"

class FluxJumpIndicator : public InternalSideIndicator
{
public:
  FluxJumpIndicator(const std::string & name, InputParameters parameters);
  virtual ~FluxJumpIndicator(){};

protected:

  virtual Real computeQpIntegral();
  void computeIndicator();
  void finalize();


};

template<>
InputParameters validParams<FluxJumpIndicator>();
