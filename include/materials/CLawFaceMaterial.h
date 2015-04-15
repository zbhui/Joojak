
#pragma once

#include "Material.h"
#include "CLawInterface.h"

using std::vector;

class CLawFaceMaterial :
public Material,
public CLawInterface
{
	struct QpValue
	{
		Real ul[10], ur[10], u_bar[10], u_diff[10];
		RealGradient dul[10], dur[10], duh[10];
		Real inv_flux[10], vis_flux[10], flux[10];
		RealVectorValue lift[10];
		Point normal;
		inline void disturbLeftValue(int component, Real ds)
		{
			ul[component] += ds;
		}
		inline void disturbRightValue(int component, Real ds)
		{
			ul[component] += ds;
		}

		inline void disturbLeftGradValue(int component, int beta, Real ds)
		{
			dul[component](beta) += ds;
		}

		inline void disturbRightGradValue(int component, int beta, Real ds)
		{
			dur[component](beta) += ds;
		}
	};

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
	Real _epsilon;

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

	void fillQpValue(QpValue &qp_value);
	void computeQpLift(QpValue &qp_value);
	void computeQpLift(RealVectorValue *lift, QpValue &qp_value);
	void computeQpFlux(Real *flux, QpValue &qp_value);
	void computeQpValue(RealVectorValue *flux_term, QpValue &qp_value);
	void fluxTerm(Real *flux, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	void liftOperator(RealVectorValue *lift, Real *ul, Real *ur);
	void addPenalty();

};

template<>
InputParameters validParams<CLawFaceMaterial>();
