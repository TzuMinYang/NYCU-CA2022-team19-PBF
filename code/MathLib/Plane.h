#pragma once

#include <MathLib/V3D.h>
#include <MathLib/P3D.h>

class Plane{
public:
	//a plane is defined by its normal, and a point that lies on it
	V3D n;
	P3D p;
public:
	Plane(void);
	Plane(const P3D& p, const V3D& n);
	Plane(const P3D& p1, const P3D& p2, const P3D& p3);
	~Plane(void);

	double getSignedDistanceToPoint(const P3D &pt) const;
	P3D getProjectionOf(const P3D &p) const;

	Plane& operator = (const Plane& other);

	//get the coefficients that define the cartesian equation of the plane: ax + by + cz + d = 0
	void getCartesianEquationCoefficients(double& a, double& b, double& c, double& d);

};
