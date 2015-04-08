
#include "CLawProblem.h"

#include "MooseApp.h"

template<>
InputParameters validParams<CLawProblem>()
{
  InputParameters params = validParams<FEProblem>();
  return params;
}

CLawProblem::CLawProblem(const std::string & name, InputParameters params) :
    FEProblem(name, params)
{
	std::cout << "配置文件：" << _app.getInputFileName() << std::endl;
	std::cout << params <<std::endl;

//    _mesh.printInfo();
}

int CLawProblem::equationIndex(const std::string& var_name)
{
	mooseError("CLawProblem::equationIndex不可调用，需要子类填充.");
}

