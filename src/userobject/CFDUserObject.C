
#include "CFDUserObject.h"

template<>
InputParameters validParams<CFDUserObject>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addClassDescription("Holds information on the cfd variable names");
  params.addRequiredCoupledVar("cfd_vars", "List of variables that represent the porepressures or (porepressure, saturations).  In single-phase models you will just have one (eg \'pressure\'), in two-phase models you will have two (eg \'p_water p_gas\', or \'p_water s_water\', etc.  These names must also be used in your kernels and material.");

  return params;
}

CFDUserObject::CFDUserObject(const std::string & name, InputParameters parameters) :
    GeneralUserObject(name, parameters),
    Coupleable(parameters, false),
    ZeroInterface(parameters),
    _num_var(coupledComponents("cfd_vars")),
    _the_names(std::string())
{
  unsigned int max_moose_var_num_seen = 0;

  _moose_var_index.resize(_num_var);
  _moose_var_value.resize(_num_var);
  _moose_var_value_old.resize(_num_var);
  _moose_grad_var.resize(_num_var);
  _moose_raw_var.resize(_num_var);
  for (unsigned int i = 0; i < _num_var; ++i)
  {
    _moose_var_index[i] = coupled("cfd_vars", i);
    max_moose_var_num_seen = (max_moose_var_num_seen > _moose_var_index[i] ? max_moose_var_num_seen : _moose_var_index[i]);
    _moose_var_value[i] = &coupledValue("cfd_vars", i); // coupledValue returns a reference (an alias) to a VariableValue, and the & turns it into a pointer
    _moose_var_value_old[i] = (_is_transient ? &coupledValueOld("cfd_vars", i) : &_zero);
    _moose_grad_var[i] = &coupledGradient("cfd_vars", i);
    _moose_raw_var[i] = getVar("cfd_vars", i);
    _the_names += getVar("cfd_vars", i)->name() + " ";
  }
  _the_names.erase(_the_names.end() - 1, _the_names.end()); // remove trailing space

  _cfd_var_index.resize(max_moose_var_num_seen + 1);
  for (unsigned int i = 0 ; i < max_moose_var_num_seen + 1 ; ++i)
	  _cfd_var_index[i] = _num_var; // NOTE: indicates that i is not a richards variable
  for (unsigned int i=0 ; i<_num_var; ++i)
	  _cfd_var_index[_moose_var_index[i]] = i;
}


unsigned int CFDUserObject::numberCFDVariables() const
{
  return _num_var;
}

unsigned int CFDUserObject::mooseVariableIndex(unsigned int cfd_var_num) const
{
  if (cfd_var_num >= _moose_var_index.size())
    mooseError("The cfd variable number " << cfd_var_num << " is out of bounds according to the RichardsVarNames UserObject");

  return _moose_var_index[cfd_var_num];
}

unsigned int CFDUserObject::cfdVariableIndex(unsigned int moose_var_num) const
{
  if (moose_var_num >= _cfd_var_index.size() || _cfd_var_index[moose_var_num] == _num_var)
    mooseError("The moose variable with number " << moose_var_num << " is not a richards according to the RichardsVarNames UserObject");

  return _cfd_var_index[moose_var_num];
}

bool CFDUserObject::notCFDVariable(unsigned int moose_var_num) const
{
  if (moose_var_num >= _cfd_var_index.size() || _cfd_var_index[moose_var_num] == _num_var)
    return true;
  return false;
}

std::string CFDUserObject::cfdVariablesNames() const
{
  return _the_names;
}

VariableValue * CFDUserObject::value(unsigned int cfd_var_num) const
{
  return _moose_var_value[cfd_var_num]; // moose_var_value is a vector of pointers to VariableValuees
}

std::vector<VariableValue *>  CFDUserObject::value() const
{
	return _moose_var_value;
}
VariableValue * CFDUserObject::valueOld(unsigned int cfd_var_num) const
{
  return _moose_var_value_old[cfd_var_num];
}

VariableGradient* CFDUserObject::gradValue(unsigned int cfd_var_num) const
{
  return _moose_grad_var[cfd_var_num];
}

MooseVariable* CFDUserObject::rawVariable(unsigned int cfd_var_num) const
{
  return _moose_raw_var[cfd_var_num];
}


