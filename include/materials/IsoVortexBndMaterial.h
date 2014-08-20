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

#include "EulerBndMaterial.h"
#include "IsoVortexBase.h"

class IsoVortexBndMaterial;

template<>
InputParameters validParams<IsoVortexBndMaterial>();

class IsoVortexBndMaterial :
public EulerBndMaterial,
public IsoVortexBase
{
public:
	IsoVortexBndMaterial(const std::string & name, InputParameters parameters);

protected:
	virtual void computeQpRightValue(Real *ur);
};
