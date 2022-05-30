#include "GreedyRandomizedMinimizer.h"

GreedyRandomizedMinimizer::GreedyRandomizedMinimizer(int p_maxIterations, int p_numTrialsBeforeAbort, bool p_printOutput){
	maxIterations = p_maxIterations;
	numTrialsBeforeAbort = p_numTrialsBeforeAbort;
	printOutput = p_printOutput;
}

GreedyRandomizedMinimizer::~GreedyRandomizedMinimizer(){
}

bool GreedyRandomizedMinimizer::minimize(ObjectiveFunction *function, dVector &p, double p_searchWindow, double & functionValue){
	dVector windowArray((int)p.size());
	windowArray.setConstant(p_searchWindow);
	return minimize(function, p, windowArray, 1, functionValue);
}



bool GreedyRandomizedMinimizer::minimize(ObjectiveFunction *function, dVector &p, const dVector& p_searchWindow, double windowScale, double & functionValue){
	if(printOutput){
		Logger::logPrint("Starting Greedy Randomized Search...\n");
	}

	//number of parameters...
	int N = (int) p.size();
	dVector pTest(N);

	function->setCurrentBestSolution(p);
	functionValue = function->computeValue(p);

	bool converged = false;
	for(int i=0; i < maxIterations; i++) {

		if (printOutput)
			Logger::logPrint("Starting iteration %d. Initial function value: %10.10lf\n", i, functionValue);

		bool foundABetterP = false;
		double testScore = 0;
		for (int j=0;j<numTrialsBeforeAbort;j++){
			for (int k = 0; k<N; k++){
				pTest[k] = p[k];
				//20% chance that we modify each parameter...
				if (getRandomNumberInRange(0, 1) > 0.2)
					pTest[k] += getRandomNumberInRange(-p_searchWindow[k]*windowScale, p_searchWindow[k]*windowScale);
			}
			testScore = function->computeValue(pTest);

			Logger::logPrint("\tsub step: %d -> objective value: %lf\n", j, testScore);
			if (testScore < functionValue){
				foundABetterP = true;
				break;
			}
		}
		if (foundABetterP == false){
			converged = true;
			if (printOutput)
				Logger::logPrint("Could not find any better solutions... aborting!\n");
		}else{
			p = pTest;
			functionValue = testScore;
			if (printOutput)
				Logger::logPrint("Found a better solution: %10.10lf\n", functionValue);
		}
	}

	function->setCurrentBestSolution(p);

	if(printOutput)
		Logger::logPrint("Done optimization step...\n");

	return converged;
}

