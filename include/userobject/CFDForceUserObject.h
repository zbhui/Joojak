
#pragma once

#include "SideUserObject.h"
#include "Eigen/Dense"

class CFDProblem;

using Eigen::Vector3d;
using std::vector;

class CFDForceUserObject : public SideUserObject
{
private:
	class Dynam
	{
	public:
		Dynam()
		{
			setZero();
		}

		void setZero()
		{
			form_force.zero();
			friction_force.zero();
			total_force.zero();
			momentum.zero();
		}
		Vector3d transToEigenVector(RealVectorValue &vector)
		{
		   return Vector3d(vector(0), vector(1), vector(2));
		}
		RealVectorValue transToMooseVector(Vector3d &vector)
		{
		   return RealVectorValue(vector(0), vector(1), vector(2));
		}

		void operator += (const Dynam &d)
		{
			this->form_force += d.form_force;
			this->friction_force += d.friction_force;
			this->total_force += d.total_force;
//			return *this;
		}

	public:
		RealVectorValue form_force;
		RealVectorValue friction_force;
		RealVectorValue total_force;
		RealVectorValue momentum;
		Point ref_point;
	};

public:
	CFDForceUserObject(const std::string & name, InputParameters parameters);

	virtual void initialize();
	virtual void finalize();
	virtual void execute();
	virtual void threadJoin(const UserObject & uo);

	virtual void computeIntegral(Dynam &local_dynam);


protected:
	void computeQpValue(Real *uh, RealGradient *duh);

private:


	CFDProblem &_cfd_problem;
	NonlinearSystem &_nl;
	THREAD_ID _tid;
	vector<VariableName> _variables;
	int _n_equations;
	int _qp;

	Real _ref_area;
	RealVectorValue _ref_point;

	Dynam _dynam;
	vector<VariableValue*> _uh;
	vector<VariableGradient*> _grad_uh;

	bool _momentum;
	bool _force;

};


template<>
InputParameters validParams<CFDForceUserObject>();
