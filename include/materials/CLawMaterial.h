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

class CLawMaterial :
public Material
{
public:
	CLawMaterial(const std::string & name, InputParameters parameters);

protected:
};

template<>
InputParameters validParams<CLawMaterial>();
