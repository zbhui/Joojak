#include "LinearElasticityMaterial.h"

template<>
InputParameters validParams<LinearElasticityMaterial>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("youngs_modulus", "杨氏模量");
  params.addRequiredParam<Real>("poissons_ratio", "泊松比");
  params.addRequiredCoupledVar("disp_x", "x位移");
  params.addRequiredCoupledVar("disp_y", "x位移");
  params.addCoupledVar("disp_z", 0, "x位移");
  return params;
}

LinearElasticityMaterial::LinearElasticityMaterial(const std::string & name, InputParameters parameters)
  :Material(name, parameters),
   _youngs_modulus(getParam<Real>("youngs_modulus")),
   _poissons_ratio(getParam<Real>("poissons_ratio")),
   _stress(declareProperty<RealTensor>("stress"))
{
  _grad_disp.resize(3);
  _grad_disp[0] = &coupledGradient("disp_x");
  _grad_disp[1] = &coupledGradient("disp_y");
  _grad_disp[2] = &coupledGradient("disp_z");
}

void LinearElasticityMaterial::computeProperties()
{
  Real lamda = _youngs_modulus*_poissons_ratio/(1+_poissons_ratio)/(1-2*_poissons_ratio);
  Real mu = _youngs_modulus/2./(1+_poissons_ratio);
  for (unsigned int qp=0; qp<_qrule->n_points(); qp++)
  {
	  RealTensor disp_tensor((*_grad_disp[0])[qp], (*_grad_disp[1])[qp], (*_grad_disp[2])[qp]);
	  RealTensor strain((disp_tensor+disp_tensor.transpose())/2.);
	  for (int alpha = 0; alpha < 3; ++alpha)
	  {
		  for (int beta = 0; beta < 3; ++beta)
		  {
			  _stress[qp](alpha, beta) = 2*mu*(strain(alpha, beta));
			  if(alpha == beta)
				  _stress[qp](alpha, beta) += lamda*strain.tr();
		  }
	  }
  }
}
