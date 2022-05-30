#pragma once

#include "MathLib/Matrix.h"
#include <Utils/Utils.h>

class ObjectiveFunction{
public:
	ObjectiveFunction();
	virtual ~ObjectiveFunction();

	//this should always return the current value of the objective function
	virtual double computeValue(const dVector& p) = 0;
	virtual void addHessianEntriesTo(DynamicArray<MTriplet>& hessianEntries, const dVector& p) { addEstimatedHessianEntriesTo(hessianEntries, p); }
	virtual void addGradientTo(dVector& grad, const dVector& p) { addEstimatedGradientTo(grad, p); }
	//this method gets called whenever a new best solution to the objective function is found
	virtual void setCurrentBestSolution(const dVector& p) { bestSolutionYet = p; }

	void addEstimatedGradientTo(dVector& gradient, const dVector& p);
	void addEstimatedHessianEntriesTo(DynamicArray<MTriplet>& hessianEntries, const dVector& p);

	void testGradientWithFD(const dVector& p);
	void testHessianWithFD(const dVector& p);

private:
public:
	dVector bestSolutionYet;

	double weight = 1.0;
	std::string description;
};

