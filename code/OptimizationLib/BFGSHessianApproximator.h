#pragma once

#include <cmath>
#include <vector>
#include <deque>

#include <MathLib/Matrix.h>

/**
	Implementation of Limited Memory BFGS. It approximates the hessian of a function f(x). It implements an efficient method to
	compute Hinv * v.
*/
class BFGSHessianApproximator {
public:
	BFGSHessianApproximator(int historyLength = 5) {
		this->historyLength = historyLength;
		reset();
	}

	void add_x_and_dfdx_to_history(const dVector& x, const dVector& dfdx) {
		if (x.size() != last_x.size() || dfdx.size() != last_dfdx.size()) {
			reset();
			last_x = x;
			last_dfdx = dfdx;
			return;
		}
		xDiffArray.push_back(x - last_x); last_x = x;
		dfdxDiffArray.push_back(dfdx - last_dfdx); last_dfdx = dfdx;

		if ((int)xDiffArray.size() > historyLength) {
			xDiffArray.pop_front();
			dfdxDiffArray.pop_front();
		}
	}

	/*
		L-BFGS inverse Hessian multiplication. Given a direction v (typically the gradient), this method
		computes the transformation H^{-1}*v, where H is the approximated Hessian matrix. Note that the multiplication
		happens in place, v will be overwritten.
		Based on description from https://en.wikipedia.org/wiki/Limited-memory_BFGS
	*/
	void compute_Hinv_v_inplace(dVector& v);
	dVector compute_Hinv_v(const dVector& v);

	void reset() {
		resize(last_x, 0);
		resize(last_dfdx, 0);
		xDiffArray.clear();
		dfdxDiffArray.clear();
	}

private:
	std::deque<dVector> xDiffArray;
	std::deque<dVector> dfdxDiffArray;
	dVector last_x;
	dVector last_dfdx;

	int historyLength = 10;
};

