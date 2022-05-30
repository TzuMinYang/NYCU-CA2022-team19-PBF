#pragma once

#include "ObjectiveFunction.h"

class GreedyRandomizedMinimizer{
public:
	GreedyRandomizedMinimizer(int p_maxIterations = 100, int p_numTrialsBeforeAbort = 300, bool p_printOutput = false);
	virtual ~GreedyRandomizedMinimizer();

	/**
		use gradient descent to optimize a function. Since the gradient of a function gives the direction
		of steepest descent, all one needs to do is go in that direction, but we'll use a line search to know how far to go...
	*/
	virtual bool minimize(ObjectiveFunction *function, dVector &p, double p_searchWindow, double & functionValue);
	virtual bool minimize(ObjectiveFunction *function, dVector &p, const dVector& p_searchWindow, double windowScale, double & functionValue);


	bool printOutput;
protected:
	int maxIterations;
	int numTrialsBeforeAbort;

};

