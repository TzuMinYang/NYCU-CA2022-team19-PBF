#pragma once

#include "MathLib/P3D.h"
#include "MathLib/V3D.h"

class CollisionPlane {
private:
	P3D pointOnPlane;
	V3D normal;
public:
	CollisionPlane(P3D p, V3D n);
	P3D handleCollision(P3D point);
};