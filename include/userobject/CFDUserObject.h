
#pragma once

#include "GeneralUserObject.h"
#include "Coupleable.h"
#include "ZeroInterface.h"

class CFDUserObject:
  public GeneralUserObject,
  public Coupleable,
  public ZeroInterface
{
 public:
	CFDUserObject(const std::string & name, InputParameters parameters);

  void initialize(){};
  void execute(){};
  void finalize(){};

  /// the number of porepressure variables
  unsigned int numberCFDVariables() const;

  /**
   * the MOOSE variable number
   * @param richards_var_num the richards variable number
   * eg if richards_vars = 'pwater pgas', and the variables in
   * the simulation are 'temperature pwater pgas displacement'
   * then moose_var_num(0) = 1
   */
  unsigned int mooseVariableIndex(unsigned int cfd_var_index) const;

  /**
   * the richards variable number
   * @param moose_var_num the MOOSE variable number
   * eg if richards_vars = 'pwater pgas', and the variables in
   * the simulation are 'temperature pwater pgas displacement'
   * then richards_var_num(2) = 1
   */
  unsigned int cfdVariableIndex(unsigned int moose_var_num) const;

  /**
   * returns true if moose_var_num is not a richards var
   * @param moose_var_num the MOOSE variable number
   * eg if richards_vars = 'pwater pgas', and the variables in
   * the simulation are 'temperature pwater pgas displacement'
   * then not_pressure_var(0) = true, no_pressure_var(1) = false
   */
  bool notCFDVariable(unsigned int moose_var_num) const;

  /**
   * a space-separated string of richards variable names
   * eg richards_names() = 'pwater pgas'
   */
  std::string cfdVariablesNames() const;

  /**
   * a vector of pointers to VariableValues
   * @param richards_var_num the pressure variable number
   * eg if richards_vars = 'pwater pgas', then
   * (*richards_vals(1))[qp] = pgas evaluated at quadpoint qp
   * Also richards_vals(i) = &coupledValue
   */
  VariableValue * value(unsigned int richards_var_num) const;
  std::vector<VariableValue *>  value() const;
  /**
   * a vector of pointers to old VariableValues
   * @param richards_var_num the richards variable number
   * eg if richards_vars = 'pwater pgas', then
   * (*richards_vals_old(1))[qp] = old pgas evaluated at quadpoint qp
   * Also richards_vals_old(i) = &coupledValueOld
   */
  VariableValue * valueOld(unsigned int richards_var_num) const;

  /**
   * a vector of pointers to grad(Variable)
   * @param richards_var_num the richards variable number
   * eg if richards_vars = 'pwater pgas', then
   * (*grad_var(1))[qp] = grad(pgas) evaluated at quadpoint qp
   * Also grad_var(i) = &coupledGradient
   */
  VariableGradient * gradValue(unsigned int richards_var_num) const;

  /**
   * The moose variable for the given richards_var_num
   * This is got using the getVar function.  It allows
   * direct extraction of nodal variable values
   * used in mass lumping.
   * @param richards_var_num the richards variable number
   */
  MooseVariable * rawVariable(unsigned int richards_var_num) const;

  /**
   * The nodal variable values for the given richards_var_num
   * To extract a the value of pressure variable "pvar", at
   * node i, use (*RichardsVarNames.nodal_var(pvar))[i]
   * @param richards_var_num the richards variable number
   */
//  VariableValue * nodal_var(unsigned int richards_var_num) const;

  /**
   * The old nodal variable values for the given richards_var_num
   * @param richards_var_num the richards variable number
   */
//  VariableValue * nodal_var_old(unsigned int richards_var_num) const;

  /// return the _var_types string
//  std::string var_types() const;


 protected:

  /// number of richards variables
  unsigned int _num_var;

  /// space-separated string of names of porepressure variables
  std::string _the_names;

  /// physical meaning of the variables.  Eg 'pppp' means 'all variables are pressure variables'
//  MooseEnum _var_types;

  /// _moose_var_num[i] = the moose variable number corresponding to richards variable i
  std::vector<unsigned int> _moose_var_index;

  /// _pressure_var_num[i] = the richards variable corresponding to moose variable i
  std::vector<unsigned int> _cfd_var_index;

  /// moose_var_value[i] = values of richards variable i
  std::vector<VariableValue *> _moose_var_value; // this is a vector of pointers to VariableValues

  /// moose_var_value_old[i] = old values of richards variable i
  std::vector<VariableValue *> _moose_var_value_old;

  /// moose_var_value[i] = values of richards variable i
//  std::vector<VariableValue *> _moose_nodal_var_value; // this is a vector of pointers to VariableValues

  /// moose_var_value_old[i] = old values of richards variable i
//  std::vector<VariableValue *> _moose_nodal_var_value_old;

  /// moose_grad_var[i] = gradient values of richards variable i
  std::vector<VariableGradient *> _moose_grad_var;

  /// _moose_raw_var[i] = getVar of richards variable i
  std::vector<MooseVariable *> _moose_raw_var;

};

template<>
InputParameters validParams<CFDUserObject>();
