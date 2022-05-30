#pragma once

//This header file contains useful constants and macros.
#include <vector>

#define DynamicArray std::vector
typedef unsigned int uint;

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdio.h>



/**
	The epsilon value is used for all kinds of numerical computations. For instance, when checking to see if two points are equal, we will check to
	see if they are equal, within epsilon.
*/
#define EPSILON 1e-10

#ifndef INFINITY
#define INFINITY DBL_MAX
#endif

#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

#define TINY 0.000000001

#define PI 3.14159265

/**
	This macro checks to see if the value of x is zero within epsilon: -epsilon<x<epsilon
*/
#define IS_ZERO(x) (fabs(x)<EPSILON)

#define IS_EQUAL(x, y) IS_ZERO(((x)-(y)))

/**
	Computes the value of x in radians
*/
#define RAD(x) (((x) * PI)/180.0)

/**
	And this computes the value of x in degrees
*/
#define DEG(x) ((x * 180)/PI)

#define SQR(x) ((x)*(x))

#define SGN(x) (((x)<0)?(-1):(1))

inline void boundToRange(double* v, double min, double max){
	if (*v < min) *v = min;
	if (*v > max) *v = max;
}

inline double safeACOS(double val){
	if (val<-1)
		return PI;
	if (val>1)
		return 0;
	return acos(val);
}

inline double safeASIN(double val){
	boundToRange(&val, -1, 1);
	return asin(val);
}

inline double getRandomNumberInRange(double min, double max){
	double range = max - min;
	int largeIntValue = RAND_MAX;
	int randVal = rand();
	double val = (randVal % largeIntValue) / (double)(largeIntValue-1);
	return min + val * range;
}

/*
	Draw a number from a gaussian distribution.
	To get a number with a certain mean and variance, take
	r = mean + sigma*getRandomGaussian()
*/
inline double getRandomGaussian(){
	double x1, x2, rquad;
	do {
		x1 = 2.0*getRandomNumberInRange(0,1) - 1.0;
		x2 = 2.0*getRandomNumberInRange(0,1) - 1.0;
		rquad = x1*x1 + x2*x2;
	} while (rquad >= 1 || rquad <= 0);

	double fac = sqrt(-2*log(rquad)/rquad);

	return fac*x2;
}

/*
	if v < min, this method returns 0. If v > max, it returns 1. For everything else it returns an interpolated value;
*/
inline double mapTo01Range(double v, double min, double max){
	double t = v;
	if (fabs(min - max) < TINY) return 1;
	boundToRange(&t, min, max); 
	t = (t-min)/(max-min);
	return t;
}

inline double linearlyInterpolate(double v1, double v2, double t1, double t2, double t){
	if (v1 == v2)
		return v2;
	return (t-t1)/(t2-t1) * v2 + (t2-t)/(t2-t1) * v1;
}

inline int roundToInt(double r){
    return (int)((r > 0.0) ? (r + 0.5) : (r - 0.5)); 
}

inline double randNumberIn01Range() {
	return ((double)rand() / ((double)RAND_MAX + 1));
}

inline bool isNaN(double x) { return (x != x); }

