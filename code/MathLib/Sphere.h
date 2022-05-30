#pragma once

/*
  A sphere
*/
#include "MathLib.h"
#include "P3D.h"

class SimpleSphere{
public:
  SimpleSphere(P3D center_, double radius_) : center(center_), radius(radius_) {}

  SimpleSphere(double x, double y, double z, double radius_) : center(P3D(x,y,z)), radius(radius_) {}


private:
  P3D center;
  double radius;
};
