#pragma once

#include "MathLib/Matrix.h"

class V3D;
/**
 * Point in 3d.
 */
class P3D : public Vector3d {
public:

/**
	Contructors and destructor
*/
	P3D();

	explicit P3D(double val);

	P3D(double x, double y, double z);

	P3D(double x, double y);

	explicit P3D(const V3D& v);

	P3D(const P3D& p);

	~P3D();

/**
	operators
*/

	P3D& operator = (const P3D& other);

	P3D& operator = (const V3D& other);

	bool operator == (const P3D& p) const;

	bool operator != (const P3D& p) const;
    
    bool operator < (const P3D& p) const;

	double at(int i) const;

	double& at(int i);

	//return *this + v
	P3D operator + (const V3D &v) const;

	//return *this - p
	V3D operator - (const P3D &p) const;

	//*this += v
	P3D& operator += (const V3D &v);

	P3D & operator-=(const V3D & v);

	//*this += p
	P3D& operator += (const P3D& p);

	//return *this + p
	P3D operator + (const P3D& p) const;

	//return *this / v
	P3D operator / (double val) const;

	//return *this * v
	P3D operator * (double val) const;

	// * this /= v
	P3D& operator /= (double val);

	// *this /= v
	P3D& operator *= (double val);

	// return - *this
	P3D operator - ();

/**
	useful methods
*/
	double getComponentAlong(const V3D& other);

	double getComponentAlong(const V3D& other) const;

	void setComponentAlong(const V3D& other, double val);

	void addOffsetToComponentAlong(const V3D& other, double offset);

	void scaleComponentAlong(const V3D& other, double scale);

	void boundComponentAlong(const V3D& other, double min, double max);

	// makes it all zero
	void zero();
};

// compute surface area of a triangle
double computeTriangleSurfaceArea(P3D & p0, P3D & p1, P3D & p2);


