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

#include "EulerCellMaterial.h"

template<>
InputParameters validParams<EulerCellMaterial>()
{
  InputParameters params = validParams<Material>();
  params += validParams<EulerBase>();
  params.addCoupledVar("disp_x", " disp_x");
  params.addCoupledVar("disp_y", " disp_y");
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

EulerCellMaterial::EulerCellMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		EulerBase(name, parameters),
		_invis_term(declareProperty<std::vector<RealVectorValue> >("cell_material")),
		_jacobi(declareProperty<std::vector<std::vector<RealVectorValue> > >("cell_jacobi_variable")),
		_disp_x(coupledValue("disp_x")),
        _disp_y(coupledValue("disp_y")),
		_disp_old_x(coupledValueOld("variables", 0)),
        _disp_old_y(coupledValueOld("disp_y"))
{
	_n_equations = coupledComponents("variables");
	for (int eq = 0; eq < _n_equations; ++eq)
		_uh.push_back(&coupledValue("variables", eq));
}

void EulerCellMaterial::computeQpProperties()
{
	if(_bnd) return;

	Real uh[10];
	computeQpValue(uh);
	_invis_term[_qp].resize(_n_equations);
	fluxTerm(&_invis_term[_qp][0], uh);

	RealVectorValue invis_term_new[10];
	_jacobi[_qp].resize(_n_equations);
	for (int p = 0; p < _n_equations; ++p)
		_jacobi[_qp][p].resize(_n_equations);

	for (int q = 0; q < _n_equations; ++q)
	{
		uh[q] += _ds;
		fluxTerm(invis_term_new, uh);
		for (int p = 0; p < _n_equations; ++p)
			_jacobi[_qp][p][q] = (invis_term_new[p] - _invis_term[_qp][p])/_ds;

		uh[q] -= _ds;
	}
}

void EulerCellMaterial::computeQpValue(Real* uh)
{
	for (int eq = 0; eq < _n_equations; ++eq)
		uh[eq] = (*_uh[eq])[_qp];
}

void EulerCellMaterial::fluxTerm(RealVectorValue* flux_term, Real* uh)
{
	RealVectorValue invis_term[10];
	inviscousTerm(invis_term, uh);
	RealVectorValue vel_grid(1, 1, 0);
	Real dt =_fe_problem.dt();
//	_console << _disp_old_x[4] <<std::endl;
//	_console << (_disp_x[4] - _disp_old_x[4])/dt <<std::endl;
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux_term[eq] = invis_term[eq] - uh[eq]*vel_grid;
	}
}
