
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
		_material_data(declareProperty<CLawCellMaterialData>("cell_material_data"))
{
}

void CLawCellMaterial::computeProperties()
{
	if(_bnd || _neighbor) return ;
	_claw_problem.computeCellMaterial(*this);
}

