/** ***********************************************************
*  @file
*  @brief
*  @author	刘  伟
*
*  Program:   Joojak
*  Copyright (c) 刘伟，张来平，2014，空气动力学国家重点实验室(SKLA)
*  All rights reserved.
*  ************************************************************
**/

#pragma once

#include "Material.h"
#include "MaterialProperty.h"

class EulerMaterial;

template<>
InputParameters validParams<EulerMaterial>();

/**
 * Euler流体的材料属性
 */
class EulerMaterial : public Material
{
public:
	EulerMaterial(const std::string & name, InputParameters parameters);

protected:
	Real pressure(Real *uh);
	Real enthalpy(Real *uh);
	Real temperature(Real  *uh);
	Real mach_local(Real *uh);
	Real acous(Real *uh);
	Real physicalViscosity(Real *uh);

  virtual void computeQpProperties();

  int _n_equations;
  Real _gamma;
  Real _prandtl;
  Real _reynolds;
  Real _mach;

  /// 积分点上的变量值
  std::vector<VariableValue*> _uh;

  MaterialProperty<RealVectorValue*> & _inviscous_term;

  void computeQpValue(Real *uh);
};
