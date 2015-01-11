#include "ElasticityKernel.h"

template<>
InputParameters validParams<ElasticityKernel>()
{
  InputParameters params = validParams<Kernel>();
  params.set<bool>("use_displaced_mesh") = true;

  return params;
}

ElasticityKernel::ElasticityKernel(const std::string & name, InputParameters parameters)
  :Kernel(name, parameters),
   _stress(getMaterialProperty<RealTensorValue>("stress"))
{}

Real ElasticityKernel::computeQpResidual()
{
	RealVectorValue temp(_stress[_qp](0,0), _stress[_qp](0,1), _stress[_qp](0,2));
	return temp *_grad_test[_i][_qp];
}

Real ElasticityKernel::computeQpJacobian()
{
    return 0;
}

Real ElasticityKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
    return 0;
}
