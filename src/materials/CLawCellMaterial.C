
#include "CLawCellMaterial.h"

template<>
InputParameters validParams<CLawCellMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addParam<Real>("ds", 1.490116119384766e-08, "微扰量");
  return params;
}

CLawCellMaterial::CLawCellMaterial(const std::string & name, InputParameters parameters):
		Material(name, parameters),
		CLawInterface(parameters),
		_ds(getParam<Real>("ds")),
		_flux_term(declareProperty<std::vector<RealVectorValue> >("flux_term")),
		_flux_jacobi_variable(declareProperty<std::vector<std::vector<RealVectorValue> > >("flux_term_jacobi_variable")),
		_flux_jacobi_grad_variable(declareProperty<std::vector<std::vector<RealTensorValue> > >("flux_term_jacobi_grad_variable")),
		_cell_material_data(declareProperty<CLawCellMaterialData>("cell_material_data"))
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		MooseVariable &val = getVariable(eq);
		mooseAssert(val.order() == _var_order, "变量的阶不同");

		_uh.push_back(_is_implicit ? &val.sln() : &val.slnOld());
		_grad_uh.push_back(_is_implicit ? &val.gradSln(): &val.gradSlnOld());
	}
}

void CLawCellMaterial::computeQpProperties()
{
	if(_bnd || _neighbor) return;

	resizeQpProperty();

	QpValue qp_value;
	fillQpValue(qp_value);
	computeQpValue(&_flux_term[_qp][0], qp_value);
	computeQpValue(_cell_material_data[_qp]._flux_term, qp_value);

	for (int q = 0; q < _n_equations; ++q)
	{
		qp_value.disturbValue(q, _ds);
		computeQpValue(qp_value);
		for (int p = 0; p < _n_equations; ++p)
			_flux_jacobi_variable[_qp][p][q] = (qp_value.flux_term[p] - _flux_term[_qp][p])/_ds;

		qp_value.disturbValue(q, -_ds);

		for (int beta = 0; beta < 3; ++beta)
		for (int q = 0; q < _n_equations; ++q)
		{
			qp_value.disturbGradValue(q, beta, _ds);

			computeQpValue(qp_value);
			for (int alpha = 0; alpha< 3; ++alpha)
			for (int p = 0; p < _n_equations; ++p)
			{
				_flux_jacobi_grad_variable[_qp][p][q](alpha, beta) = (qp_value.flux_term[p](alpha) - _flux_term[_qp][p](alpha))/_ds;
			}

			qp_value.disturbGradValue(q, beta, -_ds);
		}
	}
}

void CLawCellMaterial::fillQpValue(QpValue &qp_value)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		qp_value.uh[eq] = (*_uh[eq])[_qp];
		qp_value.duh[eq] = (*_grad_uh[eq])[_qp];
	}
}

void CLawCellMaterial::computeQpValue(QpValue &qp_value)
{
	QpValue invis, viscous;
	convertionTerm(invis.flux_term, qp_value.uh);
	diffusionTerm(viscous.flux_term, qp_value.uh, qp_value.duh);
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		qp_value.flux_term[eq] = invis.flux_term[eq] - viscous.flux_term[eq];
	}
}

void CLawCellMaterial::computeQpValue(RealVectorValue* flux_term, QpValue& qp_value)
{
	QpValue invis, viscous;
	convertionTerm(invis.flux_term, qp_value.uh);
	diffusionTerm(viscous.flux_term, qp_value.uh, qp_value.duh);
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux_term[eq] = invis.flux_term[eq] - viscous.flux_term[eq];
	}
}

void CLawCellMaterial::resizeQpProperty()
{
	_flux_term[_qp].resize(_n_equations);
	_flux_jacobi_variable[_qp].resize(_n_equations);
	_flux_jacobi_grad_variable[_qp].resize(_n_equations);
	for (int p = 0; p < _n_equations; ++p)
	{
		_flux_jacobi_variable[_qp][p].resize(_n_equations);
		_flux_jacobi_grad_variable[_qp][p].resize(_n_equations);
	}

}
