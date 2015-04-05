
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
}
