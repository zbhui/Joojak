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

#include "NSBndMaterial.h"
#include "CouetteFlowBase.h"

class CouetteFlowBndMaterial;

template<>
InputParameters validParams<CouetteFlowBndMaterial>();

class CouetteFlowBndMaterial :
public NSBndMaterial,
public CouetteFlowBase
{
public:
	CouetteFlowBndMaterial(const std::string & name, InputParameters parameters);

protected:
	virtual void computeQpRightValue(Real *ur,RealGradient *dur, Real *ul, RealGradient *dul);
};
