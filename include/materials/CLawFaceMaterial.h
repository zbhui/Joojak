
#pragma once

#include "Material.h"
#include "CLawInterface.h"

using std::vector;

class CLawFaceMaterial :
public Material,
public CLawInterface
{
public:
	CLawFaceMaterial(const std::string & name, InputParameters parameters);

protected:
	virtual void resizeQpProperty();
	virtual void computeQpProperties();

	const Real & _current_elem_volume;
	const Real & _neighbor_elem_volume;
	const Real & _current_side_volume;
	Real _ds;
	Real _sigma;
	int _n_equations;
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

	virtual void computeQpLeftValue(Real *ul);
	virtual void computeQpRightValue(Real *ur);
	virtual void computeQpLeftGradValue(RealGradient *dul);
	virtual void computeQpRightGradValue(RealGradient *dur);

	void fluxTerm(Real *flux, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	void liftOperator(RealVectorValue *lift, Real *ul, Real *ur);
	void addPenalty();

};

template<>
InputParameters validParams<CLawFaceMaterial>();
