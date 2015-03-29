
#pragma once

#include "MeshModifier.h"
#include "libmesh/fe.h"


class BuildSideSetFromBlock : public MeshModifier
{
public:
  BuildSideSetFromBlock(const std::string & name, InputParameters parameters);

  virtual ~BuildSideSetFromBlock(){};

  virtual void modify();

protected:
  std::string _boundary_name;
};

template<>
InputParameters validParams<BuildSideSetFromBlock>();

