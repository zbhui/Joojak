
#pragma once

#include "CLawBoundaryMaterial.h"
#include "CouetteFlowBase.h"

class CouetteFlowBndMaterial :
public CLawBoundaryMaterial,
public CouetteFlowBase
{
public:
	CouetteFlowBndMaterial(const std::string & name, InputParameters parameters);

protected:
	void computeQpRightValue(Real *ur);
	virtual void computeQpFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul);

};

template<>
InputParameters validParams<CouetteFlowBndMaterial>();
