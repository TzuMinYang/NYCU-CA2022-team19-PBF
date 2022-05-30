#include "Quaternion.h"
#include <Utils/Utils.h>

//	A constructor. Intuitive...
Quaternion::Quaternion(double s, const V3D& v) {
	this->s = s;
	this->v = v;
}

// A copy constructor
Quaternion::Quaternion(const Quaternion& other) {
	this->s = other.s;
	this->v = other.v;
}

// Default constructor - no rotation
Quaternion::Quaternion() {
	this->s = 1;
	this->v = V3D(0, 0, 0);
}

// Another constructor.
Quaternion::Quaternion(double w, double x, double y, double z) {
	this->s = w;
	this->v = V3D(x, y, z);
}

// destructor.
Quaternion::~Quaternion(void) {
}

// A copy operator
Quaternion& Quaternion::operator = (const Quaternion &rhs) {
	this->s = rhs.s;
	this->v = rhs.v;
	return *this;
}

// returns the result of multiplying the current quaternion by rhs
Quaternion Quaternion::operator * (const Quaternion &rhs) const{
	Quaternion q = *this;
	q *= rhs;
	return q;
}

// multiply the current quaternion by the rhs one.
Quaternion& Quaternion::operator *= (const Quaternion &rhs) {
	double newS = this->s * rhs.s - this->v.dot(rhs.v);
	V3D newV = rhs.v * this->s + this->v * rhs.s + this->v.cross(rhs.v);
	this->s = newS;
	this->v = newV;
	return *this;
}

// multiplies the quaternion by a scalar
Quaternion& Quaternion::operator *= (double scalar) {
	this->s *= scalar;
	this->v *= scalar;
	return *this;
}

// divides the quaternion by a scalar
Quaternion& Quaternion::operator /= (double scalar) {
	this->s /= scalar;
	this->v /= scalar;
	return *this;
}

// return a copy of the quaternion, multiplied by a scalar
Quaternion Quaternion::operator * (double scalar) const {
	return Quaternion(s * scalar, v * scalar);
}

// returns the result of dividing the current quaternion by rhs
Quaternion Quaternion::operator / (double scalar) const {
	return Quaternion(s / scalar, v * scalar);
}

// returns sum of *this and rhs
Quaternion Quaternion::operator + (const Quaternion &rhs) const {
	return Quaternion(s + rhs.s, v + rhs.v);
}

// adds rhs to this quaternion
Quaternion& Quaternion::operator += (const Quaternion &rhs) {
	this->s += rhs.s;
	this->v += rhs.v;
	return *this;
}

// rotates vector v by this quaternion
V3D Quaternion::operator * (const V3D &v) const {
	return rotate(v);
}

// normalize the quaternion
Quaternion& Quaternion::toUnit() {
	*this /= this->length();
	return *this;
}

// returns the rotation angle represented by this quaternion - in the range -pi to pi. Since q and -q are the same rotation, we need to know which side is "up" explicitly
double Quaternion::getRotationAngle(const V3D& positiveRotAxis) {
	int sinSign = SGN(positiveRotAxis.dot(v));
	double result = 2 * safeACOS(s);
	if (sinSign < 0)
		result = -result;
	if (result > PI) result -= 2 * PI;
	if (result < -PI) result += 2 * PI;
	return result;
}

// Returns the complex conjugate of the current quaternion.
Quaternion Quaternion::getComplexConjugate() const {
	return Quaternion(s, -v);
}

//	Returns the inverse of the current quaternion: q * q^-1 = identity quaternion: s = 1, v = (0,0,0)
Quaternion Quaternion::getInverse() const {
	return this->getComplexConjugate() / length2();
}

//	Returns the length of this quaternion.
double Quaternion::length() const {
	return sqrt(s*s + v.dot(v));
}

double Quaternion::length2() const {
	return s*s + v.dot(v);
}

// Computes the dot product between the current quaternion and the one given as parameter.
double Quaternion::dot(const Quaternion &other) const {
	return (this->s * other.s + this->v.dot(other.v));
}

// linear interpolation (don't do it!). If t = 0, the result is *this. If t = 1, result is other.
Quaternion Quaternion::linearlyInterpolateWith(const Quaternion &other, double t) const {
	if (t<0) t = 0;
	if (t>1) t = 1;
	Quaternion result = (*this)*(1 - t) + other*t;
	return result * (1 / result.length());
}

// SLERP (spherical linear interpolation). If t = 0, the result is *this. If t = 1, result is other.
Quaternion Quaternion::sphericallyInterpolateWith(const Quaternion &other, double t) const {
	Quaternion result;
	sphericallyInterpolateWith(other, t, result);
	return result;
}

// SLERP (spherical linear interpolation). If t = 0, the result is *this. If t = 1, result is other.
void Quaternion::sphericallyInterpolateWith(const Quaternion &other, double t, Quaternion &result) const {
	//make sure that we return the same value if either of the quaternions involved is q or -q 
	if (this->dot(other) < 0) {
		Quaternion temp;
		temp.s = -other.s;
		temp.v = -other.v;
		this->sphericallyInterpolateWith(temp, t, result);
		return;
	}

	if (t<0) t = 0;
	if (t>1) t = 1;
	double dotProduct = this->dot(other);
	double sinTheta = sqrt(MAX(0, 1 - dotProduct*dotProduct));
	if (sinTheta == 0)
		result = *this;
	else{
		double invSinTheta = 1.0 / sinTheta;
		double theta = safeACOS(dotProduct);
		double fN = sin(theta * (1 - t)) * invSinTheta;
		double fP = sin(theta * t) * invSinTheta;

		result.s = this->s*fN + other.s*fP;
		result.v = this->v*fN + other.v*fP;
	}
}

// return a rotated version of u
V3D Quaternion::rotate(const V3D& u) const {
	V3D res;
	rotate(u, &res);
	return res;
}

// rotate u by the inverse of the quaternion
V3D Quaternion::inverseRotate(const V3D& u) const {
	V3D res;
	inverseRotate(u, &res);
	return res;
}

// target will be set to u rotated by this quaternion. 'target' cannot point at u
void Quaternion::rotate(const V3D& u, V3D* target) const {
	//uRot = q * (0, u) * q' = (s, v) * (0, u) * (s, -v)
	//working it out we get:
	V3D t = v.cross(u) + u * s;
	*target = v.cross(t) + t*s + v*u.dot(v);
}

// target will be set to u rotated by this quaternion. 'target' cannot point at u
void Quaternion::inverseRotate(const V3D& u, V3D *target) const {
	//uRot = q * (0, u) * q' = (s, -v) * (0, u) * (s, v)
	//working it out we get:
	V3D t = u.cross(v) + u * s;
	*target = t.cross(v) + t * s + v * u.dot(v);
}

// This method returns the quaternion in axis-angle representation.
void Quaternion::getAxisAngle(V3D &axis, double &angle) const {
	if (IS_EQUAL(s, 1.0)){
		angle = 0.0;
		axis = V3D(1,0,0);
	}else{
		angle = 2 * safeACOS(s);
		axis = v * (1.0 / sqrt(1 - s*s));
		axis.toUnit();
//		assert(IS_ZERO(axis.length() - 1));
	}
}

// sets value of rotation quaternion from arguments
void Quaternion::setRotationFrom(double angle, const V3D& axis) {
	s = cos(angle / 2.0);
	v = axis * sin(angle / 2.0);
}

// sets value of rotation quaternion from arguments - the magnitude of the vector encodes the rotation angle
void Quaternion::setRotationFrom(const V3D& vec) {
	double len = vec.length();
	if (IS_ZERO(len)) {
		s = 1;
		v = V3D();
	}
	else
		setRotationFrom(len, vec / len);
}

// sets value of rotation quaternion from coordinate frame passed in as arguments. The vectors passed in are the x, y and z axis of a local 
//coordinate frame expressed in the global coordinate frame. The resulting quaternion will convert vectors from the local frame to the global frame
void Quaternion::setRotationFrom(const V3D& xRot, const V3D& yRot, const V3D& zRot) {
	assert(IS_ZERO(xRot.length()-1) && IS_ZERO(yRot.length()-1) && IS_ZERO(zRot.length()-1) && IS_ZERO(xRot.dot(yRot)) && IS_ZERO(yRot.dot(zRot)) && IS_ZERO(zRot.dot(xRot)));

	double mat[9] = {xRot[0], yRot[0], zRot[0], 
					 xRot[1], yRot[1], zRot[1],
					 xRot[2], yRot[2], zRot[2]};

#define MAT(i, j) (mat[i*3+j])

	double tr = MAT(0, 0) + MAT(1, 1) + MAT(2, 2);

	if (tr > 0) {
		double S = sqrt(tr + 1) * 2;
		this->s = 0.25*S;
		this->v[0] = (MAT(2, 1) - MAT(1, 2)) / S;
		this->v[1] = (MAT(0, 2) - MAT(2, 0)) / S;
		this->v[2] = (MAT(1, 0) - MAT(0, 1)) / S;
	}
	else if (MAT(0, 0)>MAT(1, 1) && MAT(0, 0)>MAT(2, 2)) {
		double S = sqrt(1 + MAT(0, 0) - MAT(1, 1) - MAT(2, 2)) * 2;
		this->s = (MAT(2, 1) - MAT(1, 2)) / S;
		this->v[0] = 0.25*S;
		this->v[1] = (MAT(0, 1) + MAT(1, 0)) / S;
		this->v[2] = (MAT(0, 2) + MAT(2, 0)) / S;
	}
	else if (MAT(1, 1) > MAT(2, 2)) {
		double S = sqrt(1 + MAT(1, 1) - MAT(0, 0) - MAT(2, 2)) * 2;
		this->s = (MAT(0, 2) - MAT(2, 0)) / S;
		this->v[0] = (MAT(0, 1) + MAT(1, 0)) / S;
		this->v[1] = 0.25*S;
		this->v[2] = (MAT(1, 2) + MAT(2, 1)) / S;
	}
	else {
		double S = sqrt(1 + MAT(2, 2) - MAT(0, 0) - MAT(1, 1)) * 2;
		this->s = (MAT(1, 0) - MAT(0, 1)) / S;
		this->v[0] = (MAT(0, 2) + MAT(2, 0)) / S;
		this->v[1] = (MAT(1, 2) + MAT(2, 1)) / S;
		this->v[2] = 0.25*S;
	}

	//now test the result...
	assert(this->rotate(V3D(1, 0, 0)) == xRot && this->rotate(V3D(0, 1, 0)) == yRot && this->rotate(V3D(0, 0, 1)) == zRot);
}

//sets the 3x3 top-left corner of a matrix with a rotation cooresponding to this quaternion
Matrix3x3 Quaternion::getRotationMatrix(){
	double w = s, x = v[0], y = v[1], z = v[2];

	Matrix3x3 res;
	res << 1 - 2 * y*y - 2 * z*z,	2 * x*y - 2 * w*z,	2 * x*z + 2 * w*y,
		2 * x*y + 2 * w*z,	1 - 2 * x*x - 2 * z*z,	2 * y*z - 2 * w*x,
		2 * x*z - 2 * w*y,	2 * y*z + 2 * w*x,	1 - 2 * x*x - 2 * y*y;
	return res;
}

//sets value of rotation quaternion from the matrix that is passed in as a parameter
void Quaternion::setRotationFrom(const Matrix3x3& mat){
	double tr = mat(0, 0) + mat(1, 1) + mat(2, 2);

	if (tr > 0) {
		double S = sqrt(tr + 1) * 2;
		s = 0.25*S;
		v[0] = (mat(2, 1) - mat(1, 2)) / S;
		v[1] = (mat(0, 2) - mat(2, 0)) / S;
		v[2] = (mat(1, 0) - mat(0, 1)) / S;
	}
	else if (mat(0, 0)>mat(1, 1) && mat(0, 0)>mat(2, 2)) {
		double S = sqrt(1 + mat(0, 0) - mat(1, 1) - mat(2, 2)) * 2;
		s = (mat(2, 1) - mat(1, 2)) / S;
		v[0] = 0.25*S;
		v[1] = (mat(0, 1) + mat(1, 0)) / S;
		v[2] = (mat(0, 2) + mat(2, 0)) / S;
	}
	else if (mat(1, 1) > mat(2, 2)) {
		double S = sqrt(1 + mat(1, 1) - mat(0, 0) - mat(2, 2)) * 2;
		s = (mat(0, 2) - mat(2, 0)) / S;
		v[0] = (mat(0, 1) + mat(1, 0)) / S;
		v[1] = 0.25*S;
		v[2] = (mat(1, 2) + mat(2, 1)) / S;
	}
	else {
		double S = sqrt(1 + mat(2, 2) - mat(0, 0) - mat(1, 1)) * 2;
		s = (mat(1, 0) - mat(0, 1)) / S;
		v[0] = (mat(0, 2) + mat(2, 0)) / S;
		v[1] = (mat(1, 2) + mat(2, 1)) / S;
		v[2] = 0.25*S;
	}
}

Quaternion getRotationAxisThatAlignsVectors(const V3D& a, const V3D& b) {
	V3D aUnit = a.unit();
	V3D bUnit = b.unit();
	//we first need a rotation that gets dir to be aligned with the y-axis...
	V3D rotAxis = aUnit.cross(bUnit);
	if (IS_ZERO(rotAxis.length()))
		return Quaternion();
	rotAxis.toUnit();
	double rotAngle = a.angleWith(b, rotAxis);
	return getRotationQuaternion(rotAngle, rotAxis);
}

