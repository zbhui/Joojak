
#include "CFDForceUserObject.h"

template<>
InputParameters validParams<CFDForceUserObject>()
{
  InputParameters params = validParams<SideUserObject>();
  return params;
}

CFDForceUserObject::CFDForceUserObject(const std::string & name, InputParameters parameters) :
		SideUserObject(name, parameters)
{}


void CFDForceUserObject::execute()
{
//	std::cout << "test\n";
}

void CFDForceUserObject::finalize()
{
	std::cout << "test\n";
}
