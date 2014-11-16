
#include "PressureAndSkinFrictionCoeff.h"

template<>
InputParameters validParams<PressureAndSkinFrictionCoeff>()
{
	InputParameters params = validParams<SideValueSampler>();
	params += validParams<NSBase>();
	return params;
}

PressureAndSkinFrictionCoeff::PressureAndSkinFrictionCoeff(const std::string & name, InputParameters parameters) :
			SideValueSampler(name, parameters),
			NSBase(name, parameters)
{
	_n_equations = coupledComponents("variable");
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		_uh.push_back(&coupledValue("variable", eq));
		_grad_uh.push_back(&coupledGradient("variable", eq));
	}

	std::vector<std::string> var_names;
	var_names.push_back("Cp");
	var_names.push_back("Cf");
	_values.resize(2);

	SamplerBase::setupVariables(var_names);
}

void PressureAndSkinFrictionCoeff::execute()
{
	Real uh[10];
	RealGradient duh[10];

	Matrix3d tau;
	Vector3d normal;
	Real tau_w;
	for (unsigned int _qp=0; _qp<_qrule->n_points(); _qp++)
	{
		computeQpValue(uh, duh);
		Real pre = pressure(uh);
		Real pre_inf = pressureInfity();
		_values[0] = (pre-pre_inf)/(0.5);
		//	  _values[0] = uh[0];
		for (int alpha = 0; alpha < 3; ++alpha)
			normal(alpha) = -_normals[_qp](alpha);

		stressTerm(tau, uh, duh);
		tau_w = (tau*normal).norm()/0.5;
		_values[1] = tau_w;

		SamplerBase::addSample(_q_point[_qp], _current_elem->id(), _values);

		Real vel_star = sqrt(tau_w/uh[0]);
	}
}

void PressureAndSkinFrictionCoeff::computeQpValue(Real* uh, RealGradient *duh)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		uh[eq] = (*_uh[eq])[_qp];
		duh[eq] = (*_grad_uh[eq])[_qp];
	}
}
