
#pragma once

#include "CLawMaterial.h"

using std::vector;
class CLawProblem;

class CLawFaceMaterialData
{
public:
	Real _flux[10];
	Real _flux_jacobi_variable_ee[10][10];
	Real _flux_jacobi_variable_en[10][10];
	RealGradient _flux_jacobi_grad_variable_ee;
	RealGradient _flux_jacobi_grad_variable_en;
	RealVectorValue _lift[10];
	RealVectorValue _lift_jacobi_variable[10][10];
	RealVectorValue _lift_jacobi_variable_neighbor[10][10];
};

class CLawFaceMaterial : public CLawMaterial
{
public:
	CLawFaceMaterial(const std::string & name, InputParameters parameters);

public:
	const MooseArray<Point> & normals() {return _normals;}

public:
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

	MaterialProperty<CLawFaceMaterialData> &_face_material_data;

//	virtual void computeQpProperties();
	virtual void computeProperties();
	void computeQpValue(Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	void computeQpFlux(Real *flux, RealVectorValue *lift, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
};

template<>
InputParameters validParams<CLawFaceMaterial>();
