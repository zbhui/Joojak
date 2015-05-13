
#include "CLawProblem.h"
#include "CLawCellMaterial.h"
#include "CLawFaceMaterial.h"
#include "CLawBoundaryMaterial.h"

template<>
InputParameters validParams<CLawProblem>()
{
  InputParameters params = validParams<FEProblem>();
  params.addCoupledVar("aux_variables", "耦合变量");
  return params;
}

CLawProblem::CLawProblem(const std::string & name, InputParameters params) :
    FEProblem(name, params),
	_variables(_nl.getVariableNames()),
	_n_equations(_variables.size())
{
	std::cout << "配置文件：" << _app.getInputFileName() << std::endl;
	std::cout << params <<std::endl;
	std::cout << _n_equations <<std::endl;
//    _mesh.printInfo();
}



void CLawProblem::init()
{
	FEProblem::init();
	_variables = _nl.getVariableNames();
	_n_equations = (_variables.size());
	std::cout << _n_equations <<std::endl;
}
void CLawProblem::computeFaceFlux(Real* flux, RealVectorValue* lift, Real* ul, Real* ur, RealGradient* dul, RealGradient* dur, Point& normal, Real penalty)
{
	RealVectorValue ifl[10], ifr[10], vfl[10], vfr[10];
	computeLift(lift, ul, ur, normal);

	viscousTerm(vfl, ul, dul);
	viscousTerm(vfr, ur, dur);

	fluxRiemann(flux, ul, ur, normal);

	for (int eq = 0; eq < _n_equations; ++eq)
	{
		flux[eq] -= 0.5*((vfl[eq]+vfr[eq])-penalty*lift[eq])*normal;
	}

}

void CLawProblem::computeLift(RealVectorValue *lift, Real *ul, Real *ur, Point &normal)
{
	RealGradient duh[10];
	Real uh[10];
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		duh[eq] = (ul[eq]-ur[eq])/2.*normal;
		uh[eq] = (ul[eq]+ur[eq])/2;
	}
	viscousTerm(lift, uh, duh);
}

void CLawProblem::computeBoundaryFlux(Real* flux, RealVectorValue* lift, Real* ul, RealGradient* dul, Point& normal, Real penalty, std::string bc_type)
{
	Real ur[10];
	RealGradient dur[10];
	RealVectorValue ifl[10], ifr[10], vfl[10], vfr[10];

	for (int eq = 0; eq < _n_equations; ++eq)
		dur[eq] = dul[eq];

	boundaryCondition(ur, ul, normal, bc_type);

	computeFaceFlux(flux, lift, ul, ur, dul, dur, normal, penalty);
}

void CLawProblem::computeCellFlux(RealGradient* flux, Real* uh, RealGradient* duh)
{
	RealVectorValue inv_term[10], vis_term[10], source_term[10];
	inviscousTerm(inv_term, uh);
	viscousTerm(vis_term, uh, duh);
//	sourceTerm(source_term, uh, duh);
	for (int eq = 0; eq < _n_equations; ++eq)
		flux[eq] = inv_term[eq] - vis_term[eq];
}

void CLawProblem::inviscousTerm(RealVectorValue* inviscous_term, Real* uh)
{
	mooseError("CLawProblem::inviscousTerm，需要子类填充.");

}

void CLawProblem::sourceTerm(RealVectorValue* source_term, Real* uh, RealGradient* duh)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		source_term[eq](0) = 0;
		source_term[eq](1) = 0;
		source_term[eq](2) = 0;
	}
}

void CLawProblem::viscousTerm(RealVectorValue* viscous_term, Real* uh, RealGradient* duh)
{
	mooseError("CLawProblem::viscousTerm，需要子类填充.");
}

void CLawProblem::fluxRiemann(Real* flux, Real* ul, Real* ur, Point& normal)
{
	mooseError("CLawProblem::fluxRiemann，需要子类填充.");
}

void CLawProblem::boundaryCondition(Real* ur, Real* ul, Point& normal, std::string bc_type)
{
	mooseError("CLawProblem::boundaryCondition，需要子类填充.");
}

Real CLawProblem::initialCondition(int eq)
{
	mooseError("CLawProblem::initialCondition，需要子类填充.");
}

void CLawProblem::computeCellMaterial(CLawCellMaterial& claw_cell_material)
{
	vector<VariableValue*> var = claw_cell_material._uh;
	int n_qpoints = claw_cell_material.numPoints();

	Real uh[10];
	RealVectorValue flux_term[10];
	for(int qp = 0 ; qp < n_qpoints; ++qp)
	{
		for (int eq = 0; eq < _n_equations; ++eq)
			uh[eq] =  (*var[eq])[qp];

		computeCellFlux(claw_cell_material._cell_material_data[qp]._flux_term, uh, NULL);

	}


}

void CLawProblem::computeFaceMaterial(CLawFaceMaterial& claw_face_material)
{
	vector<VariableValue*> &varl = claw_face_material._ul;
	vector<VariableValue*> &varr = claw_face_material._ur;
	int n_qpoints = claw_face_material.numPoints();

	Real _ds = 1e-08;

//	Real h_face = (_current_elem_volume+_neighbor_elem_volume)/_current_side_volume /2.;
//	Real penalty = _sigma*_var_order*_var_order/h_face;

	for(int qp = 0 ; qp < n_qpoints; ++qp)
	{
		Real ul[10], ur[10], ur_new[10];
		RealGradient dul[10], dur[10];
		Real flux_new[10];
		RealVectorValue lift_new[10];
		RealVectorValue vis_term_left[10], vis_term_right[10], vis_term_new[10];

		for (int eq = 0; eq < _n_equations; ++eq)
		{
			ul[eq] =  (*varl[eq])[qp];
			ur[eq] =  (*varr[eq])[qp];
		}

		Point normal = claw_face_material.normals()[qp];
		computeFaceFlux(claw_face_material._face_material_data[qp]._flux, claw_face_material._face_material_data[qp]._lift, ul, ur, dul, dur, normal, 0);
//		fluxRiemann(claw_face_material._face_material_data[qp]._flux, ul, ur, normal);
//		for (int q = 0; q < _n_equations; ++q)
//		{
//			ul[q] += _ds;
//			computeQpFlux(flux_new, lift_new, ul, ur, dul, dur);
//			for (int p = 0; p < _n_equations; ++p)
//			{
//				Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
//				_flux_jacobi_variable_ee[_qp][p][q] = tmp;
//				_lift_jacobi_variable[_qp][p][q] = (lift_new[p] - _lift[_qp][p])/_ds;
//			}
//			ul[q] -= _ds;
//
//			ur[q] += _ds;
//			computeQpFlux(flux_new, lift_new, ul, ur, dul, dur);
//			for (int p = 0; p < _n_equations; ++p)
//			{
//				Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
//				_flux_jacobi_variable_en[_qp][p][q] = tmp;
//				_lift_jacobi_variable_neighbor[_qp][p][q] = (lift_new[p] - _lift[_qp][p])/_ds;
//			}
//			ur[q] -= _ds;
//		computeCellFlux(claw_cell_material._cell_material_data[qp]._flux_term, uh, NULL);
	}
}
