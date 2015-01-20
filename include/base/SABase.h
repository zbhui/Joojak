
#pragma once

#include "NSBase.h"

class SABase;

template<>
InputParameters validParams<SABase>();

class SABase:
public NSBase
{
public:
	SABase(const std::string & name, InputParameters parameters);
	virtual ~SABase(){};

protected:
	int equationIndex(const std::string &var_name);
	Real eddyViscosity(Real *uh);
	virtual void inviscousTerm(RealVectorValue* inviscous_term, Real* uh);
	virtual void viscousAndSourceTerm(RealVectorValue *viscous_term, Real* source_term, Real *uh, RealGradient *duh, Real d);
	virtual void viscousTerm(RealVectorValue *viscous_term, Real *uh, RealGradient *duh, Real* source_term = NULL, Real d = 0);
//	virtual void sourceTerm(Real* source_term, Real* uh, RealGradient* duh, Real d);
protected:
	Real _cb1, _cb2, _sigma_sa, _kappa;
	Real _cw1, _cw2, _cw3, _cv1, _cv2, _cv3;
	Real _ct1, _ct2, _ct3, _ct4;
	Real _prandtl_turb;

	Real _cw3_pow6;

	Real _nu_infty;

};

