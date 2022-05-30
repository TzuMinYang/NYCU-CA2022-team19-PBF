#include "SoftUnilateralConstraint.h"

SoftUnilateralConstraint::SoftUnilateralConstraint(double stiffness, double epsilon) {
	this->epsilon = epsilon;
	a1 = stiffness;
	b1 = -0.5 * a1 * epsilon;
	c1 = -1.0 / 3 * (-b1 - a1 * epsilon) * epsilon - 1.0 / 2 * a1 * epsilon * epsilon - b1 * epsilon;

	a2 = (-b1 - a1 * epsilon) / (epsilon * epsilon);
	b2 = a1;
	c2 = b1;
	d2 = c1;
}

SoftUnilateralConstraint::~SoftUnilateralConstraint() {
}

// returns 1/2 C'C, where C is the current set of equality constraint values
double SoftUnilateralConstraint::computeValue(double x) {
	if (x < 0)
		return 0.5 * a1 * x * x + b1 * x + c1;
	if (x < epsilon)
		return 1.0 / 3 * a2 * x * x * x + 0.5 * b2 * x * x + c2 * x + d2;
	return 0;
}

double SoftUnilateralConstraint::computeDerivative(double x) {
	if (x < 0)
		return a1 * x + b1;
	if (x < epsilon)
		return a2 * x * x + b2 * x + c2;
	return 0;
}

double SoftUnilateralConstraint::computeSecondDerivative(double x) {
	if (x < 0)
		return a1;
	if (x < epsilon)
		return 2 * a2 * x + b2;
	return 0;
}


