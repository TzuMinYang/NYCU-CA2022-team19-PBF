#pragma once

#include <MathLib/P3D.h>
#include <MathLib/V3D.h>
#include <MathLib/Plane.h>

class Ray {
public:
	//a ray, or half-line, is defined by a point (where the ray starts) and a vector (direction of the ray)
	P3D origin;
	V3D direction;

public:
	Ray(void);
	Ray(const P3D &origin, const V3D &direction) : origin(origin), direction(direction.unit()) {}
	~Ray();

	// returns the point on the ray with distance 't' from the origin
	P3D getPointAt(double t) const;
	// returns the 't' value away from origin where p can be found
	double getRayParameterFor(const P3D& p) const;
	
	// returns the smallest distance between the line segment and the ray. 
	// If ClosestPtOnRay is not NULL, it will be set to the point on the ray that is closest to the segment
	double getDistanceToSegment(const P3D &p1, const P3D &p2, P3D *closestPtOnRay = NULL) const;

	// returns the distance between the line segment and the point p
	// If ClosestPtOnRay is not NULL, it will be set to the point on the ray that is closest to c
	double getDistanceToPoint(const P3D& p, P3D *closestPtOnRay = NULL) const;

	// returns the smallest distance between the line segment and the ray. 
	// and if ClosestPtOnRay is not NULL, it will be set to the point on the ray that is closest to c
	double getDistanceToPlane(const Plane& plane, P3D *closestPtOnRay = NULL) const;

	//draws the ray
	void draw();

};
