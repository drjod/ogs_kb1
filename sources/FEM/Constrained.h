
#ifndef CONSTRAINED_H_
#define CONSTRAINED_H_

#include "FEMEnums.h"

struct Constrained {
	double constrainedValue;
	FiniteElement::ProcessType constrainedProcessType;
	FiniteElement::PrimaryVariable constrainedPrimVar;
	ConstrainedType::type constrainedDirection;
	ConstrainedVariable::type constrainedVariable;
	bool _isCompleteConstrained;
	bool _completeConstrainedStateOff;
	std::vector<bool>_constrainedNodes;

	Constrained ()
	{
		constrainedValue=0.0;
		constrainedProcessType = FiniteElement::INVALID_PROCESS;
		constrainedPrimVar = FiniteElement::INVALID_PV;
		constrainedDirection = ConstrainedType::INVALID_CONSTRAINED_TYPE;
		constrainedVariable = ConstrainedVariable::INVALID_CONSTRAINED_VARIABLE;
		_isCompleteConstrained = false;
		_completeConstrainedStateOff = false;
	}

};

#endif
