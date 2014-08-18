/** ***********************************************************
*  @file
*  @brief
*  @author	刘  伟
*
*  Program:   Joojak
*  Copyright (c) 刘伟，张来平，2014，空气动力学国家重点实验室(SKLA)
*  All rights reserved.
*  ************************************************************
**/

#include "EulerBndMaterial.h"
//#include "MooseObject.h"

template<>
InputParameters validParams<EulerBndMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<CFDBase>();

  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "多个求解变量");

  return params;
}

EulerBndMaterial::EulerBndMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		CFDBase(name, parameters),
		_invis_term(declareProperty<std::vector<RealVectorValue> >("invis_term_left")),
		_invis_term_neighbor(declareProperty<std::vector<RealVectorValue> >("invis_term_right")),
		_flux_diff(declareProperty<Real>("flux_diff"))
{
	 std::vector<NonlinearVariableName> _variables = getParam<std::vector<NonlinearVariableName> >("variables");
	_n_equations = coupledComponents("variables");
	for (size_t i = 0; i < _n_equations; ++i)
	{
		SystemBase &_sys = (*parameters.get<SystemBase *>("_sys"));
		MooseVariable &val = _sys.getVariable(_tid, _variables[i]);

		_ul.push_back(_is_implicit ? &val.sln() : &val.slnOld());
//		_ur.push_back(_is_implicit ? &val.slnNeighbor() : &val.slnOldNeighbor());
	}
}

void EulerBndMaterial::computeQpProperties()
{
	Real ul[10], ur[10];
	computeQpLeftValue(ul);
	computeQpRightValue(ur);
	inviscousTerm(_invis_term[_qp], ul);
	inviscousTerm(_invis_term_neighbor[_qp], ur);
	_flux_diff[_qp] = 1.;
}

void EulerBndMaterial::computeQpLeftValue(Real* ul)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ul[eq] = (*_ul[eq])[_qp];
}

void EulerBndMaterial::computeQpRightValue(Real* ur)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		ur[eq] = (*_ul[eq])[_qp];
}
