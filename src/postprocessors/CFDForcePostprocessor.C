

#include "CFDForcePostprocessor.h"

template<>
InputParameters validParams<CFDForcePostprocessor>()
{
	MooseEnum direction_options("x y z");
	MooseEnum force_options("form friction total");
	InputParameters params = validParams<SideIntegralPostprocessor>();
	params += validParams<NSBase>();
	params.addRequiredParam<MooseEnum>("direction_by", direction_options, "气动力分量");
	params.addRequiredParam<MooseEnum>("force_type", force_options, "气动力类型");
	params.addRequiredCoupledVar("variables", "守恒变量");
	return params;
}

CFDForcePostprocessor::CFDForcePostprocessor(const std::string & name, InputParameters parameters) :
		SideIntegralPostprocessor(name, parameters),
		NSBase(name, parameters),
		_direction(getParam<MooseEnum>("direction_by")),
		_force_type(getParam<MooseEnum>("force_type"))
{
	_n_equations = coupledComponents("variables");
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		_uh.push_back(&coupledValue("variables", eq));
		_grad_uh.push_back(&coupledGradient("variables", eq));
	}
}

void CFDForcePostprocessor::threadJoin(const UserObject &y)
{
  const CFDForcePostprocessor & pps = static_cast<const CFDForcePostprocessor &>(y);
  _integral_value += pps._integral_value;
}

Real CFDForcePostprocessor::computeQpIntegral()
{
	Real uh[10];
	RealGradient duh[10];
	computeQpValue(uh, duh);

	Real pre = pressure(uh);
	Matrix3d tau;
//	Point normal = -_normals[_qp]; // 物面的外法向量和网格外向量相反

	Vector3d normal;
	stressTerm(tau, uh, duh);
	for (int alpha = 0; alpha < 3; ++alpha)
		normal(alpha) = -_normals[_qp](alpha);

	Vector3d form_force, friciton_force, total_force;
	form_force = -pre*normal/(0.5*_ref_area);
	friciton_force = tau*normal/(0.5*_ref_area);
	total_force = form_force + friciton_force;

	Quaterniond q = earthFromWind().inverse();
	form_force = q*form_force;
	friciton_force = q*friciton_force;
	total_force = q*total_force;

	switch (_force_type)
	{
	case 0: // 形状阻力
		return form_force(_direction);
	break;
	case 1: // 摩擦阻力
		return friciton_force(_direction);
	break;
	case 2: // 总阻力
		return total_force(_direction);
	break;
	default :
		return 0.;
	break;
	}

	return 0.;
}
void CFDForcePostprocessor::computeQpValue(Real* uh, RealGradient *duh)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		uh[eq] = (*_uh[eq])[_qp];
		duh[eq] = (*_grad_uh[eq])[_qp];
	}
}
