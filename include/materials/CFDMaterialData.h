
#include "libmesh/libmesh_common.h"
#include "Moose.h"

class CFDProblem;

class CFDMaterialData
{
public:
	Real r, p, t, h, s, c, q, re;
	RealVectorValue vel, mom;
	Real vel_size, vel_div, mom_size;

	RealGradient grad_rho, grad_enthalpy;
	RealTensor grad_mom, grad_vel, tau;

	Real vis;
	RealVectorValue invis_flux[10], vis_flux[10], flux[10];
	Real _gamma, _reynolds, _prandtl, _mach;
	Real uh[10];
	RealGradient duh[10];

	void reinit(CFDProblem &cfd_problem);

	virtual void InvisFlux(RealVectorValue* inviscous_term);
};
