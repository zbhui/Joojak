
#pragma once

#include "Kernel.h"
#include "Material.h"

#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"

class ElasticityKernel : public Kernel
{
public:

  ElasticityKernel(const std::string & name, InputParameters parameters);

protected:
  const MaterialProperty<RealTensorValue> & _stress;
  int _component;

  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
};

template<>
InputParameters validParams<ElasticityKernel>();
