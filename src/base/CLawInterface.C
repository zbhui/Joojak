
#include "CLawInterface.h"
#include "CLawProblem.h"

CLawInterface::CLawInterface(InputParameters& parameter):
	_claw_problem(static_cast<CLawProblem&>(*parameter.get<FEProblem *>("_fe_problem")))
{
}

void CLawInterface::convertionTerm()
{
//	_claw_base.convertionTerm();
}

void CLawInterface::diffusionTerm()
{
//	_claw_base.diffusionTerm();
}

void CLawInterface::sourceTerm()
{
//	_claw_base.sourceTerm();
}
