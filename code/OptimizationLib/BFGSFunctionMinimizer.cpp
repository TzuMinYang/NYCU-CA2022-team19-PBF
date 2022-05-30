#include "BFGSFunctionMinimizer.h"

/**
	Since the gradient of a function gives the direction of steepest descent, all one needs to do is go in that direction...
*/
void BFGSFunctionMinimizer::computeSearchDirection(ObjectiveFunction *function, const dVector &p, dVector& dp) {
	//update the hessian approximator first...
	gradient.setZero();
	function->addGradientTo(gradient, p);

	bfgsApproximator.add_x_and_dfdx_to_history(p, gradient);
	dp = bfgsApproximator.compute_Hinv_v(gradient);
}


