
#include "AddVariableAction.h"

class AddMultiVariableAction :
	public AddVariableAction
{
public:
	AddMultiVariableAction(const std::string & name, InputParameters params);
	virtual void act();

protected:
	void setInitialCondition(std::string var_name, int eq);
	void addKernel(std::string var_name, int eq);
	void addDGKernel(std::string var_name, int eq);
	void addBoundaryCondition(std::string var_name, int eq);

	std::vector<NonlinearVariableName> _variables;

};

template<>
InputParameters validParams<AddMultiVariableAction>();
