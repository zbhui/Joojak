
#include "NavierStokesProblem.h"

#include "MooseApp.h"

template<>
InputParameters validParams<NavierStokesProblem>()
{
  InputParameters params = validParams<CLawProblem>();
  params.addParam<Real>("mach",  0.2, "马赫数");
  params.addParam<Real>("gamma", 1.4, "比热比");
  params.addParam<Real>("reynolds", 1, "雷诺数");
  params.addParam<Real>("prandtl", 0.72, "prandtl数");
  params.addParam<Real>("attack", 0., "攻角");
  params.addParam<Real>("sideslip", 0., "侧滑角");
  params.addParam<Real>("pitch", 0., "俯仰角");
  params.addParam<Real>("yaw", 180., "偏航角");
  params.addParam<Real>("roll", -90., "滚转角");

  params.addParam<Real>("ref_length", 1, "参考长度");
  params.addParam<Real>("ref_area", 1, "参考面积");

  return params;
}

NavierStokesProblem::NavierStokesProblem(const std::string & name, InputParameters params) :
	CLawProblem(name, params),
	_mach(getParam<Real>("mach")),
	_gamma(getParam<Real>("gamma")),
	_reynolds(getParam<Real>("reynolds")),
	_prandtl(getParam<Real>("prandtl")),

	_attack(getParam<Real>("attack")*libMesh::pi/180),
	_sideslip(getParam<Real>("sideslip")*libMesh::pi/180),
	_pitch(getParam<Real>("pitch")*libMesh::pi/180),
	_yaw(getParam<Real>("yaw")*libMesh::pi/180),
	_roll(getParam<Real>("roll")*libMesh::pi/180),

	_ref_length(getParam<Real>("ref_length")),
	_ref_area(getParam<Real>("ref_area"))
{
}

