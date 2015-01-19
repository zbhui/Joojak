#include "ElasticityKernel.h"

template<>
InputParameters validParams<ElasticityKernel>()
{
  InputParameters params = validParams<Kernel>();
  params.set<bool>("use_displaced_mesh") = true;
  params.addRequiredParam<int>("component", "位移分量");
  return params;
}

ElasticityKernel::ElasticityKernel(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _stress(getMaterialProperty<RealTensorValue>("stress")),
   _component(getParam<int>("component"))
{}

Real ElasticityKernel::computeQpResidual()
{
	RealVectorValue stress(_stress[_qp](_component,0), _stress[_qp](_component,1), _stress[_qp](_component,2));
	return stress *_grad_test[_i][_qp];
}

Real ElasticityKernel::computeQpJacobian()
{
    return 1;
}

Real ElasticityKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
    return 0;
}
