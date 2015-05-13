
#include "CLawCellMaterial.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawCellMaterial>()
{
  InputParameters params = validParams<CLawMaterial>();
  params.addRequiredCoupledVar("variables", "守恒变量");
  return params;
}

CLawCellMaterial::CLawCellMaterial(const std::string & name, InputParameters parameter):
		CLawMaterial(name, parameter),
		_cell_material_data(declareProperty<CLawCellMaterialData>("cell_material_data"))
{
	if(_bnd || _neighbor) return ;

	for (int eq = 0; eq < _n_variables; ++eq)
	{
		_uh.push_back(&coupledValue("variables", eq));
		_grad_uh.push_back(&coupledGradient("variables", eq));
	}
}

void CLawCellMaterial::computeProperties()
{
	if(_bnd || _neighbor) return ;
	_claw_problem.computeCellMaterial(*this);
}

