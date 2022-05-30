#include "Plane.h"

Plane::Plane(void){
}

Plane::Plane(const P3D& p, const V3D& n){
	this->n = n;
	//assume the normals are unit vectors.
	this->n.toUnit();
	this->p = p;
}

Plane::Plane(const P3D& p1, const P3D& p2, const P3D& p3){
	this->p = p1;
	this->n = V3D(p1, p2).cross(V3D(p2, p3)).unit();
}

Plane::~Plane(void){
}

double Plane::getSignedDistanceToPoint(const P3D &pt) const{
	return V3D(p,pt).dot(n);
}

P3D Plane::getProjectionOf(const P3D &pt) const {
	return pt + n * -getSignedDistanceToPoint(pt);
}

Plane& Plane::operator = (const Plane& other) {
	n = other.n;
	p = other.p;
	return *this;
}

//get the coefficients that define the cartesian equation of the plane: ax + by + cz + d = 0
void Plane::getCartesianEquationCoefficients(double& a, double& b, double& c, double& d) {
	a = n[0];
	b = n[1];
	c = n[2];
	d = -(n[0]*p[0] + n[1]*p[1] + n[2]*p[2]);
}
