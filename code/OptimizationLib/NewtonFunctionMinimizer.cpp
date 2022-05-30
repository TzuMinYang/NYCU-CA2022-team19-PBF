
#include "NewtonFunctionMinimizer.h"

//TODO: try out different linear solvers, with prefactorization and without. Have a way of specifying which one this newton solver should be using...
//TODO: add suitesparse to the project as well for modified cholesky - get eigen to use it.
//TODO: get a version of the objective functions (and everything else) that works with dense matrices as well...


// The search direction is given by -Hinv * g
void NewtonFunctionMinimizer::computeSearchDirection(ObjectiveFunction *function, const dVector &p, dVector& dp) {
	timerN.restart();
	resize(H, p.size(), p.size());
	hessianEntries.clear();

	if (printOutput)
		Logger::logPrint("Time to zero out hessian: %lf\n", timerN.timeEllapsed());

	timerN.restart();
	function->addHessianEntriesTo(hessianEntries, pi);
	if (printOutput)
		Logger::logPrint("Time to compute hessian entries: %lf\n", timerN.timeEllapsed());

	timerN.restart();
	H.setFromTriplets(hessianEntries.begin(), hessianEntries.end());
	if (printOutput)
		Logger::logPrint("Time to write out hessian entries: %lf\n", timerN.timeEllapsed());

	timerN.restart();
	resize(gradient, p.size());
	function->addGradientTo(gradient, pi);
	if (printOutput)
		Logger::logPrint("Time to compute gradient: %lf\n", timerN.timeEllapsed());

	//TODO
/*
	if (svdAnalysis) {
		SVDFactorization svd;
		SparseMatrix U, S, V;
		svd.computeSVD(*function->getHessianAt(pi), U, S, V);
		for (int i = 0; i<min(S.getNumCols(), S.getNumRows()); i++)
		Logger::printStatic("%d th singluer value of hessian : %e \n", i + 1, S.getElementAt(i, i));
	}
*/

	timerN.restart();

	//dp = Hes^-1 * grad
	Eigen::SimplicialLDLT<SparseMatrix, Eigen::Lower> solver;
	//	Eigen::SparseLU<SparseMatrix> solver;
	solver.compute(H);
	dp = solver.solve(gradient);
	//dp = H.triangularView<Eigen::Lower>().solve(gradient);

	if (printOutput)
		Logger::logPrint("Time to solve linear system: %lf\n", timerN.timeEllapsed());

	int checkLinearSystemSolve = 1;
	if (checkLinearSystemSolve) {
		print("../out/hes.m", H);
		print("../out/grad.m", gradient);
		print("../out/result.m", dp);

		dVector pn = dp;
		// check how well we solved the linear system
		dVector rhs((int)gradient.size());
//		H.triangularView<Eigen::StrictlyLower>().transpose().addTo(tmpH);
		rhs = H.triangularView<Eigen::Lower>() * dp + H.transpose().triangularView<Eigen::StrictlyUpper>() * dp;
		print("../out/hes_x_dp.m", rhs);
		rhs -= gradient;
		double residual = rhs.norm();
		if (printOutput)
			Logger::logPrint("Checking linear system. Residual: %10.10lf\n", residual);
	}
}
