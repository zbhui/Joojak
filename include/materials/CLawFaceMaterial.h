
#pragma once

#include "Material.h"
#include "CLawInterface.h"

using std::vector;

class CLawProblem;

class CLawFaceMaterialData
{
public:
	void update(CLawProblem & claw_problem);
	void setProblem(CLawProblem & claw_problem, Real ds);

private:
	void computeQpValue(RealVectorValue *flux_term);
	CLawProblem * _claw_problem;
	int _n_equations;
	Real _ds;

public:
	Point normal;
	Real ul[10], ur[10], u_bar[10], u_diff[10];
	RealGradient dul[10], dur[10], duh[10];
	Real inv_flux[10], vis_flux[10], flux[10];
	Real _flux[10];
	Real _flux_jacobi_variable_ee[10][10];
	Real _flux_jacobi_variable_en[10][10];
	RealVectorValue lift[10];
	RealVectorValue _lift_jacobi_variable[10][10];
	RealVectorValue _lift_jacobi_variable_neighbor[10][10];
};

class CLawFaceMaterial :
public Material,
public CLawInterface
{
//	struct QpValue
//	{
//		Real ul[10], ur[10], u_bar[10], u_diff[10];
//		RealGradient dul[10], dur[10], duh[10];
//		Real inv_flux[10], vis_flux[10], flux[10];
//		RealVectorValue lift[10];
//		Point normal;
//		inline void disturbLeftValue(int component, Real ds)
//		{
//			ul[component] += ds;
//		}
//		inline void disturbRightValue(int component, Real ds)
//		{
//			ul[component] += ds;
//		}
//
//		inline void disturbLeftGradValue(int component, int beta, Real ds)
//		{
//			dul[component](beta) += ds;
//		}
//
//		inline void disturbRightGradValue(int component, int beta, Real ds)
//		{
//			dur[component](beta) += ds;
//		}
//	};

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
	MaterialProperty<CLawFaceMaterialData> &_face_material_data;

	void computeQpValue(Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	void computeQpFlux(Real *flux, RealVectorValue *lift, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
};

template<>
InputParameters validParams<CLawFaceMaterial>();
