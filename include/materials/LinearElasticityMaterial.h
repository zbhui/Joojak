
#pragma once

#include "Material.h"

class LinearElasticityMaterial : public Material
{
public:
  LinearElasticityMaterial(const std::string & name,
                           InputParameters parameters);

protected:
  virtual void computeProperties();

private:
  Real _youngs_modulus;
  Real _poissons_ratio;
  std::vector<VariableGradient*> _grad_disp;
  MaterialProperty<RealTensor> & _stress;
};

template<>
InputParameters validParams<LinearElasticityMaterial>();

