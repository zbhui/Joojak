
#include "CLawInterface.h"
#include "CLawBase.h"

CLawInterface::CLawInterface(InputParameters& parameter):
_claw_base(*parameter.get<CLawBase *>("_claw_base"))
{
}

void CLawInterface::convertionTerm()
{
	_claw_base.convertionTerm();
}

void CLawInterface::diffusionTerm()
{
	_claw_base.diffusionTerm();
}

void CLawInterface::sourceTerm()
{
	_claw_base.sourceTerm();
}
