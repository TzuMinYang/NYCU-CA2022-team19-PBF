#include "P3D.h"
#include "V3D.h"
#include "MathLib.h"

P3D::P3D() {
	at(0) = at(1) = at(2) = 0;
}

P3D::P3D(double val) {
	at(0) = at(1) = at(2) = val;
}

P3D::P3D(double x, double y, double z) {
	this->at(0) = x;
	this->at(1) = y;
	this->at(2) = z;
}

P3D::P3D(double x, double y) {
	this->at(0) = x;
	this->at(1) = y;
	this->at(2) = 0;
}

P3D::P3D(const V3D& v) {
	this->at(0) = v.at(0);
	this->at(1) = v.at(1);
	this->at(2) = v.at(2);
}

P3D::P3D(const P3D& p) {
	this->at(0) = p.at(0);
	this->at(1) = p.at(1);
	this->at(2) = p.at(2);
}

P3D::~P3D() {
}

P3D& P3D::operator = (const P3D& other) {
	this->at(0) = other.at(0);
	this->at(1) = other.at(1);
	this->at(2) = other.at(2);
	return *this;
}

P3D& P3D::operator = (const V3D& other) {
	this->at(0) = other.at(0);
	this->at(1) = other.at(1);
	this->at(2) = other.at(2);
	return *this;
}

bool P3D::operator == (const P3D& p) const {
	return IS_EQUAL(this->at(0), p.at(0)) && IS_EQUAL(this->at(1), p.at(1)) && IS_EQUAL(this->at(2), p.at(2));
}

bool P3D::operator != (const P3D& p) const {
	return !(*this == p);
}

bool P3D::operator < (const P3D& p) const {
    return (this->at(0) < p.at(0) - EPSILON) || (this->at(0) < p.at(0) + EPSILON && this->at(1) < p.at(1) - EPSILON)
        || (this->at(0) < p.at(0) + EPSILON && this->at(1) < p.at(1) + EPSILON && this->at(2) < p.at(2) - EPSILON);
}

double P3D::at(int i) const {
	return (*this)[i];
}

double& P3D::at(int i) {
	return (*this)[i];
}

//return *this + v
P3D P3D::operator + (const V3D &v) const {
	return P3D(this->at(0) + v.at(0), this->at(1) + v.at(1), this->at(2) + v.at(2));
}

//return *this + p
P3D P3D::operator + (const P3D& p) const{
	return P3D(this->at(0) + p.at(0), this->at(1) + p.at(1), this->at(2) + p.at(2));
}

//*this += p
P3D& P3D::operator += (const P3D& p) {
	this->at(0) += p.at(0);
	this->at(1) += p.at(1);
	this->at(2) += p.at(2);
	return *this;
}

//*this += v
P3D& P3D::operator += (const V3D &v) {
	this->at(0) += v.at(0);
	this->at(1) += v.at(1);
	this->at(2) += v.at(2);
	return *this;
}

//*this -= v
P3D& P3D::operator -= (const V3D &v) {
	this->at(0) -= v.at(0);
	this->at(1) -= v.at(1);
	this->at(2) -= v.at(2);
	return *this;
}

//return *this / v
P3D P3D::operator / (double val) const {
	return P3D(this->at(0) / val, this->at(1) / val, this->at(2) / val);
}

//return *this * v
P3D P3D::operator * (double val) const {
	return P3D(this->at(0) * val, this->at(1) * val, this->at(2) * val);
}

// * this /= v
P3D& P3D::operator /= (double val) {
	at(0) /= val;
	at(1) /= val;
	at(2) /= val;
	return *this;
}

// *this *= v
P3D& P3D::operator *= (double val) {
	at(0) *= val;
	at(1) *= val;
	at(2) *= val;
	return *this;
}

// return - *this
P3D P3D::operator - () {
	return P3D(-this->at(0), -this->at(1), -this->at(2));
}

// return the vector between the two points
V3D P3D::operator - (const P3D &p) const{
	return V3D(this->at(0) - p.at(0), this->at(1) - p.at(1), this->at(2) - p.at(2));
}

double P3D::getComponentAlong(const V3D& other) {
	return at(0) * other.at(0) + at(1) * other.at(1) + at(2) * other.at(2);
}

double P3D::getComponentAlong(const V3D& other) const {
	return at(0) * other.at(0) + at(1) * other.at(1) + at(2) * other.at(2);
}

void P3D::setComponentAlong(const V3D& other, double val) {
	double oldVal = getComponentAlong(other);
	this->at(0) += other.at(0) * (val - oldVal);
	this->at(1) += other.at(1) * (val - oldVal);
	this->at(2) += other.at(2) * (val - oldVal);
}

void P3D::addOffsetToComponentAlong(const V3D& other, double offset) {
	this->at(0) += other.at(0) * (offset);
	this->at(1) += other.at(1) * (offset);
	this->at(2) += other.at(2) * (offset);
}

void P3D::scaleComponentAlong(const V3D& other, double scale) {
	setComponentAlong(other, getComponentAlong(other) * scale);
}

void P3D::boundComponentAlong(const V3D& other, double min, double max) {
	double val = getComponentAlong(other);
	if (val < min) setComponentAlong(other, min);
	if (val > max) setComponentAlong(other, max);
}

// makes it all zero
void P3D::zero(){
	at(0) = at(1) = at(2) = 0;
}

// compute surface area of a triangle
double computeTriangleSurfaceArea(P3D & p0, P3D & p1, P3D & p2) {
	return 0.5 * V3D(p0, p1).cross(V3D(p0, p2)).length();
}



