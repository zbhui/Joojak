
#pragma once

#include "Action.h"
#include "libmesh/fe.h"

class CFDAction;

template<>
InputParameters validParams<CFDAction>();

class CFDAction : public Action
{
public:
  CFDAction(const std::string & name, InputParameters params);
  virtual void act();

protected:
  std::vector<NonlinearVariableName> _variables;
  std::vector<AuxVariableName> _aux_variables;
  std::string _init_cond_name;
  std::string _boun_cond_name;
  std::string _aux_kernel_name;
  std::string _time_kernel_name;
  std::string _cell_kernel_name;
  std::string _face_kernel_name;
  std::string _cell_material_name;
  std::string _face_material_name;

  void addVariables();
  void addAuxVariables();
  void setInitialCondition();
  void setBoundaryCondition();
  void addKernel();
  void addAuxKernel();
  void addDGKernel();
  void addMaterial();
};
