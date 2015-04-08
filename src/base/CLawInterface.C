
#include "CLawInterface.h"
#include "CLawProblem.h"

CLawInterface::CLawInterface(InputParameters& parameter):
	_claw_problem(static_cast<CLawProblem&>(*parameter.get<FEProblem *>("_fe_problem")))
{
}

void CLawInterface::convertionTerm(RealVectorValue *inviscous_term, Real *uh)
{
	_claw_problem.inviscousTerm(inviscous_term, uh);
}

void CLawInterface::diffusionTerm(RealVectorValue* viscous_term, Real* uh, RealGradient *duh)
{
	_claw_problem.viscousTerm(viscous_term, uh, duh);
}

void CLawInterface::sourceTerm()
{
//	_claw_base.sourceTerm();
}
