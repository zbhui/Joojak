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
#include "CLawInterface.h"

class CLawBoundaryMaterial :
public Material,
public CLawInterface
{
public:
	CLawBoundaryMaterial(const std::string & name, InputParameters parameters);

protected:
	std::string _bc_type;
	const Real & _current_elem_volume;
	const Real & _neighbor_elem_volume;
	const Real & _current_side_volume;
	Real _ds;
	Real _sigma;
	Real _epsilon;

	std::vector<VariableValue*> _ul;
	std::vector<VariableGradient*> _grad_ul;
	MaterialProperty<std::vector<Real> > & _flux;
	MaterialProperty<std::vector<std::vector<Real> > > &_flux_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealGradient> > > &_flux_jacobi_grad_variable;
	MaterialProperty<std::vector<RealVectorValue> > & _lift;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > > &_lift_jacobi_variable;

	virtual void computeQpProperties();
	virtual void resizeQpProperty();
	virtual void computeQpLeftValue(Real *ul, RealGradient *dul);
	void computeQpFlux(Real *flux, RealVectorValue *lift, Real *ul, RealGradient *dul);
};

template<>
InputParameters validParams<CLawBoundaryMaterial>();
