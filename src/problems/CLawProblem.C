
#include "CLawProblem.h"
#include "CLawCellMaterial.h"
#include "CLawFaceMaterial.h"
#include "CLawBoundaryMaterial.h"

template<>
InputParameters validParams<CLawProblem>()
{
  InputParameters params = validParams<FEProblem>();
  params.addParam<std::vector<VariableName> >("aux_variables", "耦合变量");
  return params;
}

CLawProblem::CLawProblem(const std::string & name, InputParameters params) :
    FEProblem(name, params),
	_aux_variables(getParam<std::vector<VariableName> >("aux_variables"))
{
	std::cout << "配置文件：" << _app.getInputFileName() << std::endl;
	std::cout << params <<std::endl;
}



void CLawProblem::init()
{
	FEProblem::init();
	_n_equations =_nl.getVariableNames().size();
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


void CLawProblem::computeCellFlux(RealGradient* flux, Real* source, Real* uh, RealGradient* duh)
{
	RealVectorValue inv_term[10], vis_term[10], source_term[10];
	inviscousTerm(inv_term, uh);
	viscousTerm(vis_term, uh, duh);
	sourceTerm(source, uh, duh);
	for (int eq = 0; eq < _n_equations; ++eq)
		flux[eq] = inv_term[eq] - vis_term[eq];
}

void CLawProblem::inviscousTerm(RealVectorValue* inviscous_term, Real* uh)
{
	mooseError("CLawProblem::inviscousTerm，需要子类填充.");

}

void CLawProblem::sourceTerm(Real* source_term, Real* uh, RealGradient* duh)
{
	for (int eq = 0; eq < _n_equations; ++eq)
	{
		source_term[eq]= 0;
		source_term[eq] = 0;
		source_term[eq] = 0;
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

Real CLawProblem::initialCondition(const Point & point, int eq)
{
	mooseError("CLawProblem::initialCondition，需要子类填充.");
}

void CLawProblem::computeCellMaterial(CLawCellMaterial& cell)
{
	vector<VariableValue*> &var = cell._uh;
	vector<VariableGradient*> &grad_var = cell._grad_uh;
	int n_qpoints = cell.numPoints();
	int n_variables = cell.numVariables();

	MaterialProperty<CLawCellMaterialData >& material = cell._material_data; ;
	Real uh[10];
	RealGradient duh[10];
	RealVectorValue flux_new[10];
	Real source_new[10];
	Real _ds = 1e-08;
	for(int qp = 0 ; qp < n_qpoints; ++qp)
	{
		for (int eq = 0; eq < n_variables; ++eq)
		{
			uh[eq] =  (*var[eq])[qp];
			duh[eq] =  (*grad_var[eq])[qp];
		}
		RealVectorValue *flux = material[qp]._flux_term;
		Real *source = material[qp]._source_term;
		computeCellFlux(flux, source, uh, duh);

		for (int q = 0; q < _n_equations; ++q)
		{
			uh[q] += _ds;
			computeCellFlux(flux_new, source_new, uh, duh);
			for (int p = 0; p < _n_equations; ++p)
			{
				material[qp]._flux_jacobi_variable[p][q] = (flux_new[p] - flux[p])/_ds;
				material[qp]._source_jacobi_variable[p][q] = (source_new[p] - source[p])/_ds;
			}
			uh[q] -= _ds;
		}

		for (int q = 0; q < _n_equations; ++q)
		for (int beta = 0; beta < 3; ++beta)
		{
			duh[q](beta) += _ds;
			computeCellFlux(flux_new, source_new, uh, duh);
			for (int p = 0; p < _n_equations; ++p)
			{
				material[qp]._source_jacobi_grad_variable[p][q](beta) = (source_new[p] - source[p])/_ds;
				for (int alpha = 0; alpha< 3; ++alpha)
				{
					material[qp]._flux_jacobi_grad_variable[p][q](alpha, beta) = (flux_new[p](alpha) - flux[p](alpha))/_ds;
				}
			}
			duh[q](beta) -= _ds;
		}
	}
}

void CLawProblem::computeFaceMaterial(CLawFaceMaterial& face)
{
	vector<VariableValue*> &varl = face._uh;
	vector<VariableValue*> &varr = face._uh_neighbor;
	vector<VariableGradient*> &grad_varl = face._grad_uh;
	vector<VariableGradient*> &grad_varr = face._grad_uh_neighbor;
	int n_qpoints = face.numPoints();
	int n_variables = face.numVariables();

	MaterialProperty<CLawFaceMaterialData >& material = face._material_data; ;

	Real _ds = 1e-08;

	for(int qp = 0 ; qp < n_qpoints; ++qp)
	{
		Real ul[10], ur[10], ur_new[10];
		RealGradient dul[10], dur[10];
		Real flux_new[10];
		RealVectorValue lift_new[10];
		RealVectorValue vis_term_left[10], vis_term_right[10], vis_term_new[10];

		for (int eq = 0; eq < n_variables; ++eq)
		{
			ul[eq] =  (*varl[eq])[qp];
			ur[eq] =  (*varr[eq])[qp];
			dul[eq] = (*grad_varl[eq])[qp];
			dur[eq] = (*grad_varr[eq])[qp];
		}

		Point normal = face.normals()[qp];
		Real penalty = face.penalty();
		Real *flux = face._material_data[qp]._flux;
		RealVectorValue *lift = face._material_data[qp]._lift;
		computeFaceFlux(flux, lift, ul, ur, dul, dur, normal, penalty);

		for (int q = 0; q < _n_equations; ++q)
		{
			ul[q] += _ds;
			computeFaceFlux(flux_new, lift_new, ul, ur, dul, dur, normal, penalty);
			for (int p = 0; p < _n_equations; ++p)
			{
				Real tmp = (flux_new[p] - flux[p])/_ds;
				material[qp]._flux_jacobi_variable_ee[p][q] = tmp;
				material[qp]._lift_jacobi_variable_ee[p][q] = (lift_new[p] - lift[p])/_ds;
			}
			ul[q] -= _ds;

			ur[q] += _ds;
			computeFaceFlux(flux_new, lift_new, ul, ur, dul, dur, normal, penalty);
			for (int p = 0; p < _n_equations; ++p)
			{
				material[qp]._flux_jacobi_variable_en[p][q] = (flux_new[p] - flux[p])/_ds;
				material[qp]._lift_jacobi_variable_en[p][q] = (lift_new[p] - lift[p])/_ds;
			}
			ur[q] -= _ds;
		}

		for (int q = 0; q < _n_equations; ++q)
		for (int beta = 0; beta < 3; ++beta)
		{
			dul[q](beta) += _ds;
			computeFaceFlux(flux_new, lift_new, ul, ur, dul, dur, normal, penalty);
			for (int p = 0; p < _n_equations; ++p)
				material[qp]._flux_jacobi_grad_variable_ee[p][q](beta) = (flux_new[p] - flux[p])/_ds;

			dul[q](beta) -= _ds;

			dur[q](beta) += _ds;
			computeFaceFlux(flux_new, lift_new, ul, ur, dul, dur, normal, penalty);
			for (int p = 0; p < _n_equations; ++p)
				material[qp]._flux_jacobi_grad_variable_en[p][q](beta) = (flux_new[p] - flux[p])/_ds;

			dur[q](beta) -= _ds;
		}
	}
}

void CLawProblem::computeBoundaryMaterial(CLawBoundaryMaterial& bnd)
{
	MaterialProperty<CLawBoundaryMaterialData>& material = bnd._material_data; ;
	vector<VariableValue*> &varl = bnd._uh;
	vector<VariableGradient*> &grad_varl = bnd._grad_uh;
	int n_qpoints = bnd.numPoints();
	int n_variables = bnd.numVariables();

	Real _ds = 1e-08;
	std::string bc_type = bnd.getBCType();

	for(_qp = 0 ; _qp < n_qpoints; ++_qp)
	{
		Real ul[10], ur[10], ur_new[10];
		RealGradient dul[10], dur[10];
		Real flux_new[10];
		RealVectorValue lift_new[10];
		RealVectorValue vis_term_left[10], vis_term_right[10], vis_term_new[10];

		for (int eq = 0; eq < n_variables; ++eq)
		{
			ul[eq] =  (*varl[eq])[_qp];
			dul[eq] = (*grad_varl[eq])[_qp];
		}

		Point normal = bnd.normals()[_qp];
		Real penalty = bnd.penalty();
		Real *flux = bnd._material_data[_qp]._flux;
		RealVectorValue *lift = bnd._material_data[_qp]._lift;

		computeBoundaryFlux(flux, lift, ul, dul, bnd);
		for (int q = 0; q < _n_equations; ++q)
		{
			ul[q] += _ds;
			computeBoundaryFlux(flux_new, lift_new, ul, dul, bnd);
			for (int p = 0; p < _n_equations; ++p)
			{
				material[_qp]._flux_jacobi_variable[p][q] = (flux_new[p] - flux[p])/_ds;
				material[_qp]._lift_jacobi_variable[p][q] = (lift_new[p] - lift[p])/_ds;
			}
			ul[q] -= _ds;
		}

		for (int q = 0; q < _n_equations; ++q)
		for (int beta = 0; beta < 3; ++beta)
		{
			dul[q](beta) += _ds;
			computeBoundaryFlux(flux_new, lift_new, ul, dul, bnd);
			for (int p = 0; p < _n_equations; ++p)
			{
				material[_qp]._flux_jacobi_grad_variable[p][q](beta) = (flux_new[p] - flux[p])/_ds;
			}
			dul[q](beta) -= _ds;
		}
	}
}

void CLawProblem::computeBoundaryFlux(Real* flux, RealVectorValue* lift, Real* ul, RealGradient* dul, CLawBoundaryMaterial& bnd)
{
	Point normal = bnd.normals()[_qp];
	Real penalty = bnd.penalty();
	RealGradient dur[10];
	Real ur[10];
	for (int eq = 0; eq < _n_equations; ++eq)
		dur[eq] = dul[eq];

	std::string bc_type = bnd.getBCType();
	boundaryCondition(ur, ul, normal, bc_type);

	computeFaceFlux(flux, lift, ul, ur, dul, dur, normal, penalty);
}

