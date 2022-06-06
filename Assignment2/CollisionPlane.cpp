#include "CollisionPlane.h"
#include "Constants.h" // add by s


CollisionPlane::CollisionPlane(P3D p, V3D n) {
	pointOnPlane = p;
	normal = n.normalized();
}

// If the given point is colliding with this plane, returns
// the projection of that point onto this plane.
// Otherwise, returns the same point.
P3D CollisionPlane::handleCollision(P3D point)
{

	// TODO: implement collision handling with planes.
	// add by s
	float distance = (point - this->pointOnPlane).dot(this->normal);
	if (distance < KERNEL_H)
		point -= this->normal * (distance - KERNEL_H); // not sure: 要不要加回KERNEL_H?

	return point;
	// end of add
}
