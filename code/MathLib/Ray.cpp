#include "Ray.h"
#include <Utils/Utils.h>

#ifdef _WIN32
#include <include/glew.h>
#else
#include <GL/glew.h>
#endif

#include <GLFW/glfw3.h>

Ray::Ray(void){
	direction = V3D(1, 0, 0);
}

Ray::~Ray(void){
}

//draws the ray
void Ray::draw() {
	glPointSize(5.0);
	glBegin(GL_POINTS);
	glVertex3d(origin[0], origin[1], origin[2]);
	glEnd();
	glPointSize(1.0);
	glBegin(GL_LINES);
	glVertex3d(origin[0], origin[1], origin[2]);
	glVertex3d(origin[0] + 100*direction[0], origin[1] + 100 * direction[1], origin[2] + 100 * direction[2]);
	glEnd();
}

// returns the point on the ray with distance 't' from the origin
P3D Ray::getPointAt(double t) const { 
	return origin + direction * t; 
}

// returns the 't' value away from origin where p can be found
double Ray::getRayParameterFor(const P3D& p) const {
	return direction.dot(p - origin);
}

// returns the smallest distance between the line segment and the ray. 
// If ClosestPtOnRay is not NULL, it will be set to the point on the ray that is closest to the segment
double Ray::getDistanceToSegment(const P3D &p1, const P3D &p2, P3D *closestPtOnRay) const {
	// (c) Dan Sunday, http://geomalgorithms.com/a07-_distance.html
	V3D u = origin + direction * 10000 - p1;
	V3D v = p2 - p1;
	V3D w = origin - p1;

	double a = u.length2();
	double b = u.dot(v);
	double c = v.length2();
	double d = u.dot(w);
	double e = v.dot(w);
	double denom = a*c - b*b;

	double sN, sD = denom;
	double tN, tD = denom;

	if (denom <= TINY) {
		// lines are almost parallel
		sN = 0;
		sD = 1;
		tN = e;
		tD = c;
	} else {
		sN = (b*e - c*d);
		tN = (a*e - b*d);
		if (sN < 0){
			// the closest point is behind the ray's starting point
			sN = 0;
			tN = e;
			tD = c;
		}
	}

	if (tN < 0) {
		// the closest point is behind t1
		tN = 0;
		if (-d < 0)
			sN = 0;
		else if (-d > a)
			sN = sD;
		else {
			sN = -d;
			sD = a;
		}
	}
	else if (tN > tD) {
		// the closest point is behind t2
		tN = tD;
		if ((-d + b) < 0)
			sN = 0;
		else if ((-d + b) > a)
			sN = sD;
		else {
			sN = -d + b;
			sD = a;
		}
	}

	double sc = (abs(sN) < TINY) ? 0 : sN / sD;
	double tc = (abs(tN) < TINY) ? 0 : tN / tD;

	if (closestPtOnRay != NULL)
		*closestPtOnRay = origin + direction * (sc * 10000);

	return (w + u*sc - v*tc).length();
}

// returns the point on the ray that is closest to c.
// If ClosestPtOnRay is not NULL, it will be set to the point on the ray that is closest to c
double Ray::getDistanceToPoint(const P3D& p, P3D *closestPtOnRay) const {
	double d = direction.dot(p - origin);
	P3D res = origin;
	if (d > 0)
		res = getPointAt(d);

	if (closestPtOnRay)
		*closestPtOnRay = res;

	return V3D(p, res).length();
}

// returns the smallest distance between the line segment and the ray. 
// and if ClosestPtOnRay is not NULL, it will be set to the point on the ray that is closest to c
double Ray::getDistanceToPlane(const Plane& plane, P3D *closestPtOnRay) const {
	//check to see if the ray is parallel to the plane...
	if (fabs(direction.dot(plane.n)) < TINY) {
		if (closestPtOnRay)
			*closestPtOnRay = origin;
		return fabs(plane.getSignedDistanceToPoint(origin));
	}

	//we know that p = origin + t * direction lies on the plane, which means that dot product between plane p and p, dotted with n is 0...
	V3D tmpV = plane.p - origin;
	double t = tmpV.dot(plane.n) / direction.dot(plane.n);

	//check the solution now...
	assert(IS_ZERO(V3D(plane.p, origin + direction * t).dot(plane.n)));

	if (t < 0) {
		if (closestPtOnRay)
			*closestPtOnRay = origin;
		return fabs(plane.getSignedDistanceToPoint(origin));
	}

	if (closestPtOnRay)
		*closestPtOnRay = origin + direction * t;

	return 0;
}

