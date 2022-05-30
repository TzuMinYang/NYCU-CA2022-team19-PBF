#pragma once
#ifdef _WIN32
#include <include/glew.h>
#else
#include <GL/glew.h>
#endif

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include <vector>
#include <string>
#include <MathLib/MathLib.h>
#include <MathLib/P3D.h>
#include <MathLib/V3D.h>
#include <MathLib/Quaternion.h>

// this method will take a screenshot of the current scene and it will save it to a file with the given name
void saveScreenShot(char* fileName, int x, int y, int width, int height);

// prints a string using openGL
void glprint(double x, double y, const char *fmt, ...);
void glprint(double x, double y, const std::string& line);
void glprint(double x, double y, const std::vector<std::string>& lines);

void drawCircle(double cx, double cy, double r, int num_segments = 10);

//drawTorus
void drawTorus(P3D center, V3D norm, double Radius, double TubeRadius, int Sides = 20, int Rings = 30);

void drawPointTrajectory(const DynamicArray<P3D>& pointTraj, const V3D& color = V3D(1,0,0), int thickness = 1.0);

// draws a sphere of radius r, centered at origin. Discretized using nPoints parameter
void drawSphere(const P3D& origin, double r, int nPoints = 6);

// draws a cylinder from origin to origin + direction
void drawCylinder(const P3D& origin, const V3D& direction, double r, int nPoints = 6, bool cappedEnds = true);

// draws a cylinder from start to end
void drawCylinder(const P3D& start, const P3D& end, double r, int nPoints = 6, bool cappedEnds = true);

// draws a capsule from origin to origin + direction
void drawCapsule(const P3D& origin, const V3D& direction, double r, int nPoints = 6);

// draws a capsule from start to end
void drawCapsule(const P3D& start, const P3D& end, double r, int nPoints = 6);

// draws a cone from org to org+v, the tip.
void drawCone(P3D org, V3D v, double r, int nrPoints = 6);

// draws an arrow from start to end
void drawArrow(const P3D& start, const P3D& end, double r, int nPoints = 6);

// draw an axis-aligned box defined by the two corner points
void drawBox(const P3D& min, const P3D& max);

//updates a menu item to include a set of items
void generateMenuEnumFromFileList(const std::string& menuName, const DynamicArray<std::string>& itemList);

void drawMOIApproximation(const Matrix3x3 MOI, double mass);

// draws a capsule from origin to origin + direction
//void drawCapsule(const P3D& origin, const V3D& direction, double r, int nPoints = 6);

// draws a capsule between the two end points
//void drawCapsule(const P3D& start, const P3D& end, double r, int nPoints = 6);


void checkOGLErrors();

