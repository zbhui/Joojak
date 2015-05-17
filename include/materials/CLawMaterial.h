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

#pragma once

#include "Material.h"
using std::vector;
class CLawProblem;

class CLawMaterial : public Material
{
public:
	CLawMaterial(const std::string & name, InputParameters parameters);

protected:
	CLawProblem &_claw_problem;
	NonlinearSystem &_nl;
	THREAD_ID _tid;
	vector<VariableName> _variables;
	vector<VariableName> _aux_variables;
	int _n_variables;
	int _var_order;

	const Real & _current_elem_volume;
	const Real & _neighbor_elem_volume;
	const Real & _current_side_volume;
public:
	const MooseArray<Point> & qpoints() {return _q_point;}
	const MooseArray<Point> & normals() {return _normals;}
	int numPoints() {return _qrule->n_points();}
	int numVariables() {return _n_variables;}

	vector<VariableValue*> _uh;
	vector<VariableValue*> _uh_neighbor;
	vector<VariableGradient*> _grad_uh;
	vector<VariableGradient*> _grad_uh_neighbor;
protected:
};

template<>
InputParameters validParams<CLawMaterial>();
