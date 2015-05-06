
#pragma once

#include "Material.h"
#include "CLawCellMaterialData.h"

using std::vector;
class CLawProblem;

class CLawCellMaterial :
public Material
{
public:
	CLawCellMaterial(const std::string & name, InputParameters parameters);

protected:
	CLawProblem &_claw_problem;
	NonlinearSystem &_nl;
	THREAD_ID _tid;
	vector<VariableName> _variables;
	vector<VariableName> _aux_variables;
	int _num_nonliner_variables;
	int _num_aux_variables;
	int _num_variables;
	int _var_order;

	Real _ds;

	std::vector<VariableValue*> _uh;
	std::vector<VariableGradient*> _grad_uh;
	MaterialProperty<CLawCellMaterialData >& _cell_material_data;

	virtual void computeQpProperties();
	void fillQpValue();
	void computeMaterialData();
};

template<>
InputParameters validParams<CLawCellMaterial>();
