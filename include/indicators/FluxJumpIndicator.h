
#pragma once

#include "InternalSideIndicator.h"
#include "TwoMaterialPropertyInterface.h"
#include "CLawFaceMaterial.h"

class FluxJumpIndicator :
public InternalSideIndicator,
protected TwoMaterialPropertyInterface
{
public:
  FluxJumpIndicator(const std::string & name, InputParameters parameters);
  virtual ~FluxJumpIndicator(){};

protected:

  virtual Real computeQpIntegral();
  void computeIndicator();
  void finalize();

  MaterialProperty<CLawFaceMaterialData> &_face;

};

template<>
InputParameters validParams<FluxJumpIndicator>();
