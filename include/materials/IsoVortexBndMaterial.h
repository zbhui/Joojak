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

#include "CLawBoundaryMaterial.h"
#include "IsoVortexBase.h"

class IsoVortexBndMaterial :
public CLawBoundaryMaterial,
public IsoVortexBase
{
public:
	IsoVortexBndMaterial(const std::string & name, InputParameters parameters);

protected:
	void computeQpRightValue(Real *ur);
	virtual void computeQpFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul);

};

template<>
InputParameters validParams<IsoVortexBndMaterial>();
