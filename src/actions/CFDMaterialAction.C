
#include "CFDMaterialAction.h"
#include "AddVariableAction.h"
#include "Factory.h"
#include "ActionFactory.h"
#include "MooseApp.h"
#include "MooseObject.h"
#include "Parser.h"
#include "FEProblem.h"

#include "libmesh/string_to_enum.h"

template<>
InputParameters validParams<CFDMaterialAction>()
{
  InputParameters params = validParams<Action>();
  params.addParam<std::vector<SubdomainName> >("block", "The list of boundary IDs from the mesh where the object will be applied");

  return params;
}

CFDMaterialAction::CFDMaterialAction(const std::string & name, InputParameters params) :
    Action(name, params)
{
}

void CFDMaterialAction::act()
{
//	InputParameters params = _factory.getValidParams("EulerCellMaterial");
//	params.set<std::vector<SubdomainName> >("block") = getParam<std::vector<SubdomainName> >("block");
//	_problem->addMaterial("EulerCellMaterial", "cell_material1", params);
}
