
#pragma once

#include "AuxKernel.h"

class NearestWallDistance;

template<>
InputParameters validParams<NearestWallDistance>();

/**
 * Coupled auxiliary value
 */
class NearestWallDistance :
public AuxKernel
{
public:
	NearestWallDistance(const std::string & name, InputParameters parameters);

protected:
    virtual Real computeValue();

protected:
    /// 势函数
    VariableValue& _psi;
    /// 势函数梯度
    VariableGradient& _grad_psi;
};

