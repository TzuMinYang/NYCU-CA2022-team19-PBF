#include "V3D.h"
#include "P3D.h"
#include <assert.h>

#include "V3D.h"
#include "P3D.h"
#include <assert.h>


V3D::V3D() {
	at(0) = at(1) = at(2) = 0;
}

V3D::V3D(double val) {
	at(0) = at(1) = at(2) = val;
}

V3D::V3D(double x, double y, double z) {
	this->at(0) = x;
	this->at(1) = y;
	this->at(2) = z;
}

V3D::V3D(double x, double y) {
	this->at(0) = x;
	this->at(1) = y;
	this->at(2) = 0;
}

V3D::V3D(const V3D& v) {
	this->at(0) = v.at(0);
	this->at(1) = v.at(1);
	this->at(2) = v.at(2);
}

V3D::V3D(const Vector3d& v) {
	this->at(0) = v[0];
	this->at(1) = v[1];
	this->at(2) = v[2];
}

V3D::V3D(const P3D& p) {
	this->at(0) = p.at(0);
	this->at(1) = p.at(1);
	this->at(2) = p.at(2);
}

V3D::V3D(const V3D &other, double scale) {
	at(0) = other.at(0) * scale;
	at(1) = other.at(1) * scale;
	at(2) = other.at(2) * scale;
}

V3D::V3D(const P3D &p1, const P3D &p2) {
	at(0) = p2.at(0) - p1.at(0);
	at(1) = p2.at(1) - p1.at(1);
	at(2) = p2.at(2) - p1.at(2);
}

V3D::~V3D() {
}

V3D& V3D::operator = (const P3D& other) {
	this->at(0) = other.at(0);
	this->at(1) = other.at(1);
	this->at(2) = other.at(2);
	return *this;
}

V3D& V3D::operator = (const V3D& other) {
	this->at(0) = other.at(0);
	this->at(1) = other.at(1);
	this->at(2) = other.at(2);
	return *this;
}

bool V3D::operator == (const V3D& v) const {
	return IS_EQUAL(this->at(0), v.at(0)) && IS_EQUAL(this->at(1), v.at(1)) && IS_EQUAL(this->at(2), v.at(2));
}

bool V3D::operator != (const V3D& v) const {
	return !(*this == v);
}

double V3D::at(int i) const {
	return (*this)[i];
}

double& V3D::at(int i) {
	return (*this)[i];
}

//return *this + v
V3D V3D::operator + (const V3D &v) const {
	return V3D(this->at(0) + v.at(0), this->at(1) + v.at(1), this->at(2) + v.at(2));
}

//*this += v
V3D& V3D::operator += (const V3D &v) {
	this->at(0) += v.at(0);
	this->at(1) += v.at(1);
	this->at(2) += v.at(2);
	return *this;
}

//*this -= v
V3D& V3D::operator -= (const V3D &v) {
	this->at(0) -= v.at(0);
	this->at(1) -= v.at(1);
	this->at(2) -= v.at(2);
	return *this;
}

//return *this / v
V3D V3D::operator / (double val) const {
	return V3D(this->at(0) / val, this->at(1) / val, this->at(2) / val);
}

//return *this * v
V3D V3D::operator * (double val) const {
	return V3D(this->at(0)*val, this->at(1)*val, this->at(2)*val);
}

// * this /= v
V3D& V3D::operator /= (double val) {
	at(0) /= val;
	at(1) /= val;
	at(2) /= val;
	return *this;
}

// *this *= v
V3D& V3D::operator *= (double val) {
	at(0) *= val;
	at(1) *= val;
	at(2) *= val;
	return *this;
}

// return - *this
V3D V3D::operator - () const {
	return V3D(-this->at(0), -this->at(1), -this->at(2));
}

// *this - v
V3D V3D::operator - (const V3D &v) const {
	return V3D(this->at(0) - v.at(0), this->at(1) - v.at(1), this->at(2) - v.at(2));
}

double V3D::getComponentAlong(const V3D& other) {
	return at(0) * other.at(0) + at(1) * other.at(1) + at(2) * other.at(2);
}

double V3D::getComponentAlong(const V3D& other) const {
	return at(0) * other.at(0) + at(1) * other.at(1) + at(2) * other.at(2);
}

void V3D::setComponentAlong(const V3D& other, double val) {
	double oldVal = getComponentAlong(other);
	this->at(0) += other.at(0) * (val - oldVal);
	this->at(1) += other.at(1) * (val - oldVal);
	this->at(2) += other.at(2) * (val - oldVal);
}

void V3D::addOffsetToComponentAlong(const V3D& other, double offset) {
	this->at(0) += other.at(0) * (offset);
	this->at(1) += other.at(1) * (offset);
	this->at(2) += other.at(2) * (offset);
}

void V3D::scaleComponentAlong(const V3D& other, double scale) {
	setComponentAlong(other, getComponentAlong(other) * scale);
}

void V3D::boundComponentAlong(const V3D& other, double min, double max) {
	double val = getComponentAlong(other);
	if (val < min) setComponentAlong(other, min);
	if (val > max) setComponentAlong(other, max);
}

// Returns a new vector obtained by rotating the current vector. Alpha is specified in radians, and axis is assumed to be a unit vector
V3D V3D::rotate(double alpha, const V3D &axis) const {
	return rotateVec(*this, alpha, axis);
}


// computes two vectors (a and b) that are orhtogonal to the current vector.
void V3D::getOrthogonalVectors(V3D& a, V3D& b) const {
	//try to choose a vector in the y-z plane, if the z-coordinate is significant enough
	if (at(0)*at(0) / (length2()) < 0.5) //if the x component of the current vector is not too large, relatively speaking, then choose x as the main component for a (they won't be aligned perfectly)
		a = V3D(1, 0, 0);
	else //if the x component of the current vector is large, then the y (or z) components cannot be the only non-zero values, so a is not conlinear with the current vector
		a = V3D(0, 1, 0);
	b = this->cross(a); //b is orthogonal to both the current vector and a
	a = b.cross(*this); //and a is orthogonal to both a and this...
	a.toUnit();
	b.toUnit();
}

// computes dot product with v
double V3D::dot(const V3D &v) const {
	return (this->at(0) * v.at(0) + this->at(1)*v.at(1) + this->at(2) * v.at(2));
}

// computes cross product with v
V3D V3D::cross(const V3D &v) const {
	/*
	Easiest way to figure it out is to set it up like this:
	___________________
	Ux | Uy   Uz   Ux   Uy | Uz
	Vx | Vy   Vz   Vx   Vy | Vz
	-------------------
	Cross product is given by cross multiplying the items in the box, and subing the other
	diagonal
	*/
	V3D result, u = *this;
	result.at(0) = u.at(1) * v.at(2) - u.at(2) * v.at(1);
	result.at(1) = u.at(2) * v.at(0) - u.at(0) * v.at(2);
	result.at(2) = u.at(0) * v.at(1) - u.at(1) * v.at(0);

	return result;
}

double sqrt__(double val) { return sqrt(val); }

// computes the length of the vector
double V3D::length() const {
	return sqrt__(length2());
}

// computes the squared length of the vector
double V3D::length2() const {
	return at(0) * at(0) + at(1) * at(1) + at(2) * at(2);
}

// computes the projection of the current vector on v
V3D V3D::getProjectionOn(const V3D &v) const {
	/*  Here's how it goes...
	proj = V/|V| * |U| * cos(angle) = V/|V| * |U| * U.V / |U| / |V|
	after simplifying we get:
	U.V
	proj = ------- * V
	|V|*|V|
	*/
	V3D u = *this;
	return v*(u.dot(v) / (v.length2()));
}

// returns the angle between this vector and v - this method only returns angles between 0 and PI. The angle returned is measured in radians.
double V3D::angleWith(const V3D &v) const {
	// U.V = |U|*|V|*cos(angle)
	// therefore angle = inverse cos (U.V/(|U|*|V|))
	double result = this->dot(v) / (this->length() * v.length());
	return safeACOS(result);
}


// returns the angle between between this vector and v. Given a direction n, the angle returned is between -PI and PI.
double V3D::angleWith(const V3D &v, const V3D& n) const {
	// U.V = |U|*|V|*cos(angle)
	double cosAVal = this->dot(v) / (this->length() * v.length());
	V3D crossProd = this->cross(v);
	// |U x V| = |U|*|V|*sin(angle)
	double sinAVal = crossProd.length() / (this->length() * v.length());

	if (crossProd.dot(n) < 0) sinAVal *= -1;

	return atan2(sinAVal, cosAVal);
}

// normalizes the vector
V3D& V3D::toUnit() {
	double d = this->length();
	if (IS_ZERO(d))
		return *this;
	*this /= d;
	return (*this);
}

// returns a unit vector
V3D V3D::unit() const {
	V3D v(this->at(0), this->at(1), this->at(2));
	v.toUnit();
	return v;
}

// makes it all zero
void V3D::zero() {
	at(0) = at(1) = at(2) = 0;
}


V3D rotateVec(const V3D& v, double alpha, const V3D &axis) {
	assert(IS_EQUAL(axis.length2(), 1));
	double xP = axis.at(0);
	double yP = axis.at(1);
	double zP = axis.at(2);
	double cosa = cos(alpha);
	double sina = sin(alpha);

	double s[3][3] = { { 0,			-zP,		yP },
	{ zP,			0,			-xP },
	{ -yP,			xP,			0 } };
	double UUT[3][3] = { { xP*xP,		xP*yP,		xP*zP },
	{ yP*xP,		yP*yP,		yP*zP },
	{ zP*xP,		zP*yP,		zP*zP } };
	double I[3][3] = { { 1,0,0 },{ 0,1,0 },{ 0,0,1 } };
	double R[3][3] = { { 0,0,0 },{ 0,0,0 },{ 0,0,0 } };

	for (int i = 0;i<3;i++)
		for (int j = 0;j<3;j++)
		R[i][j] = UUT[i][j] + cosa*(I[i][j] - UUT[i][j]) + sina*s[i][j];

	//now that we finally have the transformation matrix set up, we can rotate the vector
	V3D result;

	result.at(0) = (R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2]);
	result.at(1) = (R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2]);
	result.at(2) = (R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]);

	return result;
}

Matrix3x3 V3D::outerProductWith(const V3D &v) {
	return *((Vector3d* const)this) * v.transpose();
}