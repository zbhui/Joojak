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

#include "CLawMaterial.h"

template<>
InputParameters validParams<CLawMaterial>()
{
  InputParameters params = validParams<Material>();
  return params;
}

CLawMaterial::CLawMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		CLawInterface(parameters)
{
}
