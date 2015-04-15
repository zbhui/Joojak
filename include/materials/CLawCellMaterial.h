
#pragma once

#include "Material.h"
#include "CLawInterface.h"

class CLawCellMaterialData
{
public:
	RealVectorValue _flux_term[10];
	RealVectorValue _flux_jacobi_variable[10][10];
	RealTensorValue _flux_jacobi_grad_variable[10][10];

public:
	void print(){}
};

class CLawCellMaterial :
public Material,
public CLawInterface
{
	struct QpValue
	{
		Real uh[10];
		RealGradient duh[10];
		RealVectorValue flux_term[10];
		inline void disturbValue(int component, Real ds)
		{
			uh[component] += ds;
		}
		inline void disturbGradValue(int component, int beta, Real ds)
		{
			duh[component](beta) += ds;
		}
	};

public:
	CLawCellMaterial(const std::string & name, InputParameters parameters);

protected:
	virtual void resizeQpProperty();
	virtual void computeQpProperties();
	Real _ds;

	std::vector<VariableValue*> _uh;
	std::vector<VariableGradient*> _grad_uh;
	MaterialProperty<std::vector<RealVectorValue> > & _flux_term;
	MaterialProperty<std::vector<std::vector<RealVectorValue> > >& _flux_jacobi_variable;
	MaterialProperty<std::vector<std::vector<RealTensorValue> > >& _flux_jacobi_grad_variable;
	MaterialProperty<CLawCellMaterialData >& _cell_material_data;

	void fillQpValue(QpValue &qp_value);
	void computeQpValue(QpValue &qp_value);
	void computeQpValue(RealVectorValue *flux_term, QpValue &qp_value);
};

template<>
InputParameters validParams<CLawCellMaterial>();
