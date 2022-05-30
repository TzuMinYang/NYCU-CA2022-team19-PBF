#pragma once

#include "MathLib.h"
#include "MathLib/Matrix.h"

class P3D;
/**
* Point in 3d.
*/
class V3D : public Vector3d {
public:

/**
	Contructors and destructor
*/

	V3D();

	explicit V3D(double val);

	V3D(double x, double y, double z);

	V3D(double x, double y);

	V3D(const V3D& v);

	V3D(const Vector3d& v);

	explicit V3D(const P3D& p);

	V3D(const V3D &other, double scale);

	V3D(const P3D &p1, const P3D &p2);

	~V3D();

/**
	operators
*/

	V3D& operator = (const V3D& other);

	V3D& operator = (const P3D& other);

	bool operator == (const V3D& v) const;

	bool operator != (const V3D& v) const;

	double at(int i) const;

	double& at(int i);

	//return *this + v
	V3D operator + (const V3D &v) const;

	//return *this - v
	V3D operator - (const V3D &v) const;

	//*this += v
	V3D& operator += (const V3D &v);

	//*this -= v
	V3D& operator -= (const V3D &v);

	//return *this / v
	V3D operator / (double val) const;

	//return *this * v
	V3D operator * (double val) const;

	// * this /= v
	V3D& operator /= (double val);

	// *this /= v
	V3D& operator *= (double val);

	// return - *this
	V3D operator - () const;

/**
	useful methods
*/

	double getComponentAlong(const V3D& other);

	double getComponentAlong(const V3D& other) const;

	void setComponentAlong(const V3D& other, double val);

	void addOffsetToComponentAlong(const V3D& other, double offset);

	void scaleComponentAlong(const V3D& other, double scale);

	void boundComponentAlong(const V3D& other, double min, double max);

	// Returns a new vector obtained by rotating the current vector. Alpha is specified in radians, and axis is assumed to be a unit vector
	V3D rotate(double alpha, const V3D &axis) const;

	// computes two vectors (a and b) that are orhtogonal to the current vector.
	void getOrthogonalVectors(V3D& a, V3D& b) const;

	// computes dot product with v
	double dot(const V3D &v) const;

	// computes cross product with v
	V3D cross(const V3D &v) const;

	// computes the length of the vector
	double length() const;

	// computes the squared length of the vector
	double length2() const;

	// computes the projection of the current vector on v
	V3D getProjectionOn(const V3D &v) const;

	// returns the angle between this vector and v - this method only returns angles between 0 and PI. The angle returned is measured in radians.
	double angleWith(const V3D &v) const;

	// returns the angle between between this vector and v. Given a direction n, the angle returned is between -PI and PI.
	double angleWith(const V3D &v, const V3D& n) const;

	// normalizes the vector
	V3D& toUnit();

	// returns a unit vector
	V3D unit() const;

	// makes it all zero
	void zero();

	Matrix3x3 outerProductWith(const V3D &v);

};

// Returns a new vector obtained by rotating v. Alpha is specified in radians, and axis is assumed to be a unit vector
V3D rotateVec(const V3D& v, double alpha, const V3D &axis);

// Returns a (uniformly) random unit vector
V3D getRandomUnitVector();

