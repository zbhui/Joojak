
#include "CFDProblem.h"
#include "CFDMaterialData.h"

void CFDMaterialData::reinit(CFDProblem &cfd_problem)
{
	r = uh[0];
	p = cfd_problem.pressure(uh);
}
