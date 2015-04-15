
#pragma once

#include "Material.h"
#include "CLawInterface.h"

using std::vector;

class CLawBoundaryMaterial :
public Material,
public CLawInterface
{
public:
	CLawBoundaryMaterial(const std::string & name, InputParameters parameters);

protected:
	virtual void resizeQpProperty();
	virtual void computeQpProperties();

//	MooseEnum _bc_type;

	const Real & _current_elem_volume;
	const Real & _neighbor_elem_volume;
	const Real & _current_side_volume;
	Real _ds;
	Real _sigma;
	Real _epsilon;

	/// 积分点上的变量值
	vector<VariableValue*> _ul;
	vector<VariableValue*> _ur;
	vector<VariableGradient*> _grad_ul;
	vector<VariableGradient*> _grad_ur;

	MaterialProperty<vector<Real> > & _flux;
	MaterialProperty<vector<vector<Real> > > & _flux_jacobi_variable_ee;
	MaterialProperty<vector<vector<Real> > > & _flux_jacobi_variable_en;
	MaterialProperty<vector<vector<RealGradient> > > & _flux_jacobi_grad_variable_ee;
	MaterialProperty<vector<vector<RealGradient> > > & _flux_jacobi_grad_variable_en;

	MaterialProperty<vector<RealVectorValue> > & _lift;
	MaterialProperty<vector<vector<RealVectorValue> > > & _lift_jacobi_variable;
	MaterialProperty<vector<vector<RealVectorValue> > > & _lift_jacobi_variable_neighbor;

	void computeQpValue(Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	virtual void computeQpLeftValue(Real *ul);
	virtual void computeQpRightValue(Real *ur);
	virtual void computeQpLeftGradValue(RealGradient *dul);
	virtual void computeQpRightGradValue(RealGradient *dur);

	void fluxTerm(Real *flux, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	void liftOperator(RealVectorValue *lift, Real *ul, Real *ur);
	void addPenalty();

};

template<>
InputParameters validParams<CLawBoundaryMaterial>();
