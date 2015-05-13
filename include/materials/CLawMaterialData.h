
#pragma once
#include "MooseVariable.h"
#include "Assembly.h"
#include "InputParameters.h"
#include "NonlinearSystem.h"

using std::vector;

class CLawProblem;
class FEProblem;
class SubProblem;

class CLawMaterialData
{
public:
	CLawMaterialData(const std::string & name, InputParameters param);
	virtual ~CLawMaterialData();

public:
	SubProblem & _subproblem;
	FEProblem & _fe_problem;
	CLawProblem & _claw_problem;
	Assembly & _assembly;
	NonlinearSystem &_nl;
	vector<VariableName> _var_name;
	int _n_variables;
	vector<MooseVariable*> _moose_var;
	int _tid;
	bool _bnd;
	bool _neighbor;
//	bool _implicit;

public:
	void computeMaterial();
	bool isBoundary();
	bool isElement();
	bool isSideSet();
};

