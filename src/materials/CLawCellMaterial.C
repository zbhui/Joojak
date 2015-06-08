
#include "CLawCellMaterial.h"
#include "CLawProblem.h"

template<>
InputParameters validParams<CLawCellMaterial>()
{
  InputParameters params = validParams<CLawMaterial>();
//  params.addCoupledVar("variables", "守恒变量");
  return params;
}

CLawCellMaterial::CLawCellMaterial(const std::string & name, InputParameters parameter):
		CLawMaterial(name, parameter),
		_material_data(declareProperty<CLawCellMaterialData>("cell_material_data"))
{
//	addMooseVariableDependency(&_fe_problem.getVariable(0, "rho"));
//	addMooseVariableDependency(&_fe_problem.getVariable(0, "momentum_x"));
//	addMooseVariableDependency(&_fe_problem.getVariable(0, "momentum_y"));
//	addMooseVariableDependency(&_fe_problem.getVariable(0, "momentum_z"));
//	addMooseVariableDependency(&_fe_problem.getVariable(0, "rhoe"));
}

void CLawCellMaterial::computeProperties()
{
	if(_bnd || _neighbor) return ;
	_claw_problem.computeCellMaterial(*this);
}

