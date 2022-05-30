#pragma once

#include <math.h>
#include <MathLib/MathLib.h>
#include <MathLib/Matrix.h>
#include "ObjectiveFunction.h"

/*!
	class used to model unilateral constraints of the type x > 0 using a C2 penalty energy f(x).
*/
class SoftUnilateralConstraint {
private:
	double a1, b1, c1, a2, b2, c2, d2, epsilon;
public:
	SoftUnilateralConstraint(double stiffness, double epsilon);

	virtual ~SoftUnilateralConstraint();

	//comptue f(x)
	double computeValue(double x);

	//compute df/dx
	double computeDerivative(double x);

	//compute ddf/dxdx
	double computeSecondDerivative(double x);
};
