#pragma once

#include "MathLib/V3D.h"
#include "MathLib/P3D.h"
#include <vector>

struct Particle {
	P3D x_i;
	V3D v_i;
	P3D x_star;
	double lambda_i;
	V3D delta_p;
	double density;
	V3D vorticity_W;
	V3D vorticity_N;
	std::vector<int> neighbors;
};