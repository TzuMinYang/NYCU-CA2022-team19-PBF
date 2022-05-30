#include "GradientDescentFunctionMinimizer.h"

/**
	Since the gradient of a function gives the direction of steepest descent, all one needs to do is go in that direction...
*/
void GradientDescentFunctionMinimizer::computeSearchDirection(ObjectiveFunction *function, const dVector &p, dVector& dp) {
	dp.setZero();
	function->addGradientTo(dp, p);
}
