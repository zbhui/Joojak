
#include "CFDForceUserObject.h"
#include "CFDProblem.h"

using Eigen::Quaterniond;

template<>
InputParameters validParams<CFDForceUserObject>()
{
	InputParameters params = validParams<SideUserObject>();
	params.addParam<Real>("ref_area"  , 1, "参考面积");
	params.addParam<Point>("ref_point", "力距参考点");
	params.addParam<Real>("ref_length", "参考长度");
	return params;
}

CFDForceUserObject::CFDForceUserObject(const std::string & name, InputParameters parameters) :
	SideUserObject(name, parameters),
	_cfd_problem(static_cast<CFDProblem&>(*parameters.get<FEProblem *>("_fe_problem"))),
	_nl(_cfd_problem.getNonlinearSystem()),
	_tid(parameters.get<THREAD_ID>("_tid")),
	_variables(_nl.getVariableNames()),
	_n_equations(_variables.size()),
	_qp(0),
	_dynam(),
	_momentum(isParamValid("ref_point") ? true : false),
	_force   (isParamValid("ref_area")  ? true : false)
{
	if(_force)
		_ref_area = getParam<Real>("ref_area");

	if(_momentum)
	{
		 _ref_point = getParam<Point>("ref_point");
	}

	for (int eq = 0; eq < _n_equations; ++eq)
	{
		MooseVariable &val = _cfd_problem.getVariable(_tid, _variables[eq]);
		_uh.push_back(&val.sln());
		_grad_uh.push_back(&val.gradSln());
	}
}

void CFDForceUserObject::initialize()
{
	_dynam.setZero();  // 归零上时刻的力系数
}

void CFDForceUserObject::execute()
{
	Dynam local_dynam;
	computeIntegral(local_dynam);                   // 计算一个面上的力系数
	_dynam += local_dynam;                        // 所有面求和
}

void CFDForceUserObject::finalize()
{
	Quaterniond q = _cfd_problem._attitude.earthFromWind().inverse();
    Vector3d form_force = q*_dynam.transToEigenVector(_dynam.form_force);
    Vector3d friction_force = q*_dynam.transToEigenVector(_dynam.friction_force);
    Vector3d total_force = q*_dynam.transToEigenVector(_dynam.total_force);


	std::string line("\n***************************************************\n");
	_console << COLOR_CYAN << line << COLOR_DEFAULT;
	_console << "Form Force: " << form_force.transpose() << COLOR_DEFAULT << '\n';

	_console<< "Friction Force: " << friction_force.transpose() << '\n';

	_console<< "Total Force: " << total_force.transpose() << '\n';

	if(_momentum)
	{
		RealVectorValue mom= _dynam.total_force.cross(_ref_point);
		Vector3d momentum = q* _dynam.transToEigenVector(mom);
		_console << "Momentun: " <<  momentum.transpose() << '\n';
	}

	_console << COLOR_CYAN << line << COLOR_DEFAULT;

}

void CFDForceUserObject::threadJoin(const UserObject& uo)  // 不同分区求和
{
	const CFDForceUserObject & pps = static_cast<const CFDForceUserObject &>(uo);
	const Dynam &temp = pps._dynam;
	_dynam += pps._dynam;
}

void CFDForceUserObject::computeIntegral(Dynam &local_dynam)
{
	Real uh[10];
	RealGradient duh[10];
	for (_qp = 0; _qp < _qrule->n_points(); _qp++)
	{
		computeQpValue(uh, duh);

		RealTensor tau;
		Real pre = _cfd_problem.pressure(uh);
		_cfd_problem.stressTerm(tau, uh, duh);
		RealVectorValue normal = -_normals[_qp];

		RealVectorValue form_force, friction_force, total_force;
		form_force = -pre*normal/(0.5*_ref_area);
		friction_force = tau*normal/(0.5*_ref_area);
		total_force = form_force + friction_force;

		local_dynam.form_force += _JxW[_qp]*_coord[_qp]*form_force;
		local_dynam.friction_force += _JxW[_qp]*_coord[_qp]*friction_force;
		local_dynam.total_force += _JxW[_qp]*_coord[_qp]*total_force;
	}
}

void CFDForceUserObject::computeQpValue(Real* uh, RealGradient *duh)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		uh[eq] = (*_uh[eq])[_qp];
		duh[eq] = (*_grad_uh[eq])[_qp];
	}
}
