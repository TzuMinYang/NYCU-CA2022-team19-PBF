#include "ObjectiveFunction.h"
#include <Utils/Logger.h>

ObjectiveFunction::ObjectiveFunction(){
}

ObjectiveFunction::~ObjectiveFunction(){
}

//this is a very slow method that evaluates the hessian of the objective function through FD...
void ObjectiveFunction::addEstimatedHessianEntriesTo(DynamicArray<MTriplet>& hessianEntries, const dVector& params){
	dVector pSet = params;

	double dp = 10e-6;
	int p = (int)pSet.size();

	dVector C_P(p), C_M(p), H_i_col(p);

	for (int i=0;i<p;i++){
		C_P.setZero();
		C_M.setZero();
		double tmpVal = pSet[i];
		pSet[i] = tmpVal + dp;
		addGradientTo(C_P, pSet);

		pSet[i] = tmpVal - dp;
		addGradientTo(C_M, pSet);

		//now reset the ith param value
		pSet[i] = tmpVal;
		//and compute the row of the hessian
		H_i_col = (C_P - C_M) / (2 * dp);

		//each vector is a column vector of the hessian, so copy it in place...
		for (int j = i;j < p;j++)
			if (!IS_ZERO(H_i_col[j]))
				hessianEntries.push_back(MTriplet(j, i, H_i_col[j]));
	}
}

void ObjectiveFunction::addEstimatedGradientTo(dVector& gradient, const dVector& params){
	if (gradient.size() != params.size())
		resize(gradient, params.size());

	dVector pSet = params;

	double dp = 10e-6; //10e-6;
	int p = (int)pSet.size();

	double f_P, f_M;
	//this is a very slow method that evaluates the hessian of the objective function through FD...
	for (int i=0;i<p;i++){
		double tmpVal = pSet[i];
		pSet[i] = tmpVal + dp;
		f_P = computeValue(pSet);

		pSet[i] = tmpVal - dp;
		f_M = computeValue(pSet);

		//now reset the ith param value
		pSet[i] = tmpVal;
		gradient[i] += (f_P - f_M)/(2*dp);
	}
	//now restore the parameter set and make sure to tidy up a bit...
	computeValue(pSet);
}

void ObjectiveFunction::testGradientWithFD(const dVector& p){
	dVector FDGradient;
	dVector analyticGradient;
	addEstimatedGradientTo(FDGradient, p);
	addGradientTo(analyticGradient, p);
    Logger::logPrint("Objective Function: testing gradients...\n");
	for (int i=0;i<p.size();i++){
		double err = FDGradient[i] - analyticGradient[i];
		if (fabs(err) > 0.0001 && 2*fabs(err) / (fabs(FDGradient[i]) + fabs(analyticGradient[i])) > 0.001)
			Logger::logPrint("Mismatch element %d: Analytic val: %lf, FD val: %lf. Error: %lf\n", i, analyticGradient[i], FDGradient[i], err);
	}
}

void ObjectiveFunction::testHessianWithFD(const dVector& p){
	SparseMatrix FDHessian(p.size(), p.size());
	SparseMatrix analyticHessian(p.size(), p.size());
	DynamicArray<MTriplet> hessianEntries;

	addEstimatedHessianEntriesTo(hessianEntries, p);
	FDHessian.setFromTriplets(hessianEntries.begin(), hessianEntries.end());
	hessianEntries.clear();
	addHessianEntriesTo(hessianEntries, p);
	analyticHessian.setFromTriplets(hessianEntries.begin(), hessianEntries.end());
	Logger::logPrint("Objective Function: testing hessians...\n");
	for (int i=0;i<p.size();i++){
		for (int j=0;j<p.size();j++){
            double err = FDHessian.coeff(i, j) - analyticHessian.coeff(i, j);
			if (fabs(err) > 0.001)
				Logger::logPrint("Mismatch element %d,%d: Analytic val: %lf, FD val: %lf. Error: %lf\n", i, j, analyticHessian.coeff(i,j), FDHessian.coeff(i,j), err);
		}
	}
}






