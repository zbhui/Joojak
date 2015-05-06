
#include "AddAuxVariableAction.h"

class AddMultiAuxVariableAction :
	public AddAuxVariableAction
{
public:
	AddMultiAuxVariableAction(const std::string & name, InputParameters params);
	virtual void act();

protected:
	void addAuxVariable(std::string var_name);
	void addAuxKernel(std::string var_name);

	std::string _type;
	std::vector<AuxVariableName> _variables;
};

template<>
InputParameters validParams<AddMultiAuxVariableAction>();
