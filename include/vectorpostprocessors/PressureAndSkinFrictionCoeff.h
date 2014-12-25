
#pragma once

#include "SideValueSampler.h"
#include "NSBase.h"

class PressureAndSkinFrictionCoeff;

template<>
InputParameters validParams<PressureAndSkinFrictionCoeff>();

class PressureAndSkinFrictionCoeff :
  public SideValueSampler,
  public NSBase
{
public:
  PressureAndSkinFrictionCoeff(const std::string & name, InputParameters parameters);
  virtual ~PressureAndSkinFrictionCoeff() {}

  virtual void execute();
protected:
	int _n_equations;
	std::vector<VariableValue*> _uh;
	std::vector<VariableGradient*> _grad_uh;
	void computeQpValue(Real *uh, RealGradient *duh);
};
