
#include "CLawCellMaterial.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawCellMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addParam<Real>("ds", 1.490116119384766e-08, "微扰量");
  params.addRequiredCoupledVar("variables", "守恒变量");
  return params;
}

CLawCellMaterial::CLawCellMaterial(const std::string & name, InputParameters parameter):
		Material(name, parameter),
		_claw_problem(static_cast<CLawProblem&>(_fe_problem)),
		_nl(_claw_problem.getNonlinearSystem()),
		_tid(parameter.get<THREAD_ID>("_tid")),
		_variables(_nl.getVariableNames()),
		_n_equations(_variables.size()),
		_n_variables(coupledComponents("variables")),
		_var_order(_claw_problem.getVariable(_tid, _variables[0]).order()),

		_ds(getParam<Real>("ds")),
		_cell_material_data(declareProperty<CLawCellMaterialData>("cell_material_data"))
{
	if(_bnd || _neighbor) return ;

	for (int eq = 0; eq < _n_variables; ++eq)
	{
		_uh.push_back(&coupledValue("variables", eq));
		_grad_uh.push_back(&coupledGradient("variables", eq));
	}
}

void CLawCellMaterial::computeQpProperties()
{
	if(_bnd || _neighbor) return ;

	for (int eq = 0; eq < _n_equations; ++eq)
	{
		_cell_material_data[_qp].uh[eq] =  (*_uh[eq])[_qp];
		_cell_material_data[_qp].duh[eq] = (*_grad_uh[eq])[_qp];
		_cell_material_data[_qp].setProblem(_claw_problem, _ds);
	}

	_cell_material_data[_qp].update(_claw_problem);
}

void CLawCellMaterial::fillQpValue()
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		_cell_material_data[_qp].uh[eq] =  (*_uh[eq])[_qp];
		_cell_material_data[_qp].duh[eq] = (*_grad_uh[eq])[_qp];
	}
}

void CLawCellMaterial::computeMaterialData()
{
	_cell_material_data[_qp].update(_claw_problem);
}
