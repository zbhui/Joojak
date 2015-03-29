
#include "CLawBase.h"

CLawBase::CLawBase(const std::string& name, InputParameters& paramter) :
	_name(name),
	_subproblem(*paramter.get<SubProblem *>("_subproblem")),
	_mesh(*paramter.get<MooseMesh *>("_mesh"))
{
}
