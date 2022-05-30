#include "GLUtils.h"
#include <AntTweakBar.h>
#include <Utils/Logger.h>
#include <Utils/Image.h>
#include <Utils/BMPIO.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iosfwd>
#include <fstream>
#include "GLContentManager.h"
#include "FreeType.h"

/**
	This method will take a screenshot of the current scene and it will save it to a file with the given name
*/
void saveScreenShot(char* fileName, int x, int y, int width, int height){
	if (fileName == NULL)
		return;
	std::ofstream out(fileName, std::ofstream::binary);
	if(!out){
		Logger::print("Cannot save screenshot %s\n", fileName);
		return;
	}

	glReadBuffer(GL_BACK);

	Image *img = new Image(3, width, height, NULL);

	glReadPixels(x , y, width, height, GL_RGB, GL_UNSIGNED_BYTE, img->getDataPointer());

	BMPIO b(fileName);
	b.writeToFile(img);

	delete img;
}

void drawCircle(double cx, double cy, double r, int num_segments) {
	glBegin(GL_LINE_LOOP);
	for (int ii = 0; ii < num_segments; ii++){
		double theta = 2.0f * 3.1415926f * double(ii) / double(num_segments);//get the current angle 

		double x = r * cos(theta);//calculate the x component 
		double y = r * sin(theta);//calculate the y component 

		glVertex3d(x + cx, y + cy, 0);//output vertex 

	}
	glEnd();
}


//print in an openGL window. The raster position needs to be already defined.
void glprint(double x, double y, const char *fmt, ...) {
	char *pBuffer = NULL;
	GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);

	FreeTypeFont* font = GLContentManager::getFont("../data/fonts/arial.ttf 12");
	font->print(x, y, pBuffer);

	RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer);
}

void glprint(double x, double y, const std::string& line) {
	FreeTypeFont* font = GLContentManager::getFont("../data/fonts/arial.ttf 12");
	font->print(x, y, line);
}

void glprint(double x, double y, const std::vector<std::string>& lines) {
	FreeTypeFont* font = GLContentManager::getFont("../data/fonts/arial.ttf 12");
	font->print(x, y, lines);
}

void drawPointTrajectory(const DynamicArray<P3D>& pointTraj, const V3D& color, int thickness) {
	glLineWidth((GLfloat)thickness);
	glColor3d(color[0], color[1], color[2]);

	glBegin(GL_LINE_STRIP);
	for (uint i = 0; i < pointTraj.size(); i++)
		glVertex3d(pointTraj[i][0], pointTraj[i][1], pointTraj[i][2]);

	glEnd();
	glLineWidth(1);
}

void drawTorus(P3D center, V3D norm, double Radius, double TubeRadius, int Sides, int Rings){
	Eigen::Quaterniond q;
	Eigen::Vector3d v1, v2;
	v1 = V3D(0, 0, 1);
	v2 = norm;
	q.setFromTwoVectors(v1, v2);
	Matrix3x3 rot = q.toRotationMatrix();


	double sideDelta = 2.0 * PI / Sides;
	double ringDelta = 2.0 * PI / Rings;
	double theta = 0;
	double cosTheta = 1.0;
	double sinTheta = 0.0;

	double phi, sinPhi, cosPhi;
	double dist;

	for (int i = 0; i < Rings; i++)
	{
		double theta1 = theta + ringDelta;
		double cosTheta1 = cos(theta1);
		double sinTheta1 = sin(theta1);

		glBegin(GL_QUAD_STRIP);
		phi = 0;
		for (int j = 0; j <= Sides; j++)
		{
			phi = phi + sideDelta;
			cosPhi = cos(phi);
			sinPhi = sin(phi);
			dist = Radius + (TubeRadius * cosPhi);
			Vector3d v1, v2;
			v1 = rot * V3D(cosTheta * cosPhi, sinTheta * cosPhi, sinPhi);
			v2 = rot * V3D(cosTheta * dist, sinTheta * dist, TubeRadius * sinPhi);
			glNormal3d(v1[0], v1[1], v1[2]);
			glVertex3d(v2[0], v2[1], v2[2]);

			v1 = rot * V3D(cosTheta1 * cosPhi, sinTheta1 * cosPhi, sinPhi);
			v2 = rot * V3D(cosTheta1 * dist, sinTheta1 * dist, TubeRadius * sinPhi);
			glNormal3d(v1[0], v1[1], v1[2]);
			glVertex3d(v2[0], v2[1], v2[2]);
		}
		glEnd();
		theta = theta1;
		cosTheta = cosTheta1;
		sinTheta = sinTheta1;

	}


}

// draws a sphere of radius r, centered at origin. Discretized using nPoints parameter
void drawSphere(const P3D& origin, double r, int nPoints) {
	//this is the normal vector
	V3D n, v;

	//now we'll draw a half sphere...
	P3D p, q;
	double angle = PI / nPoints;
	p[0] = r*cos(-PI / 2);
	p[1] = r*sin(-PI / 2);
	p[2] = 0;
	for (int i = nPoints / 2;i >= 0;i--) {
		glBegin(GL_QUAD_STRIP);
		q[0] = r*cos(-i*angle);
		q[1] = r*sin(-i*angle);
		q[2] = 0;
		//make sure we compute the normal as well as the node coordinates
		n = V3D(p).toUnit();
		glNormal3d(n[0], n[1], n[2]);
		glVertex3d(origin[0] + p[0], origin[1] + p[1], origin[2] + p[2]);
		//make sure we compute the normal as well as the node coordinates
		n = V3D(q).toUnit();
		glNormal3d(n[0], n[1], n[2]);
		glVertex3d(origin[0] + q[0], origin[1] + q[1], origin[2] + q[2]);

		for (int j = 0;j <= nPoints;j++) {
			//make sure we compute the normal as well as the node coordinates
			V3D v = V3D(p[0] * cos(2 * j * angle), p[1], p[0] * sin(2 * j * angle));
			n = v.unit();
			glNormal3d(n[0], n[1], n[2]);
			glVertex3d(origin[0] + v[0], origin[1] + v[1], origin[2] + v[2]);

			//make sure we compute the normal as well as the node coordinates
			v = V3D(q[0] * cos(2 * j * angle), q[1], q[0] * sin(2 * j * angle));
			n = v.unit();
			glNormal3d(n[0], n[1], n[2]);
			glVertex3d(origin[0] + v[0], origin[1] + v[1], origin[2] + v[2]);

		}
		p = q;
		glEnd();
	}

	//and now draw the other half - horrible, quick research code...

	p[0] = r*cos(PI / 2);
	p[1] = r*sin(PI / 2);
	p[2] = 0;
	for (int i = nPoints / 2;i >= 0;i--) {
		glBegin(GL_QUAD_STRIP);
		q[0] = r*cos(i*angle);
		q[1] = r*sin(i*angle);
		q[2] = 0;
		//make sure we compute the normal as well as the node coordinates
		n = V3D(p).toUnit();
		glNormal3d(n[0], n[1], n[2]);
		glVertex3d(origin[0] + p[0], origin[1] + p[1], origin[2] + p[2]);
		//make sure we compute the normal as well as the node coordinates
		n = V3D(q).toUnit();
		glNormal3d(n[0], n[1], n[2]);
		glVertex3d(origin[0] + q[0], origin[1] + q[1], origin[2] + q[2]);

		for (int j = 0;j <= nPoints;j++) {
			//make sure we compute the normal as well as the node coordinates
			V3D v = V3D(q[0] * cos(2 * j * angle), q[1], q[0] * sin(2 * j * angle));
			n = v.unit();
			glNormal3d(n[0], n[1], n[2]);
			glVertex3d(origin[0] + v[0], origin[1] + v[1], origin[2] + v[2]);

			//make sure we compute the normal as well as the node coordinates
			v = V3D(p[0] * cos(2 * j * angle), p[1], p[0] * sin(2 * j * angle));
			n = v.unit();
			glNormal3d(n[0], n[1], n[2]);
			glVertex3d(origin[0] + v[0], origin[1] + v[1], origin[2] + v[2]);
		}
		p = q;
		glEnd();
	}
}

void drawCapsule(const P3D& origin, const V3D& direction, double r, int nPoints) {
	if (IS_ZERO(direction.length()) || IS_ZERO(r)) return;

	//we will do this in a normalized coordinate frame, since the spheres need to match the end points of the cylinder...
	P3D org = origin;
	V3D axis(0, 1, 0);
	//we first need a rotation that gets dir to be aligned with the y-axis...
	V3D rotAxis = direction.unit().cross(axis);
	double rotAngle = safeASIN(rotAxis.length());
	if (direction.dot(axis) < 0) {
		org += direction;
		rotAxis = -direction.unit().cross(axis);
		rotAngle = safeASIN(rotAxis.length());
	}
	rotAxis.toUnit();

	glPushMatrix();
	glTranslated(org[0], org[1], org[2]);
	if (rotAxis.length() > 0.9)
		glRotated(DEG(-rotAngle), rotAxis[0], rotAxis[1], rotAxis[2]);

	org = P3D();

	//we'll start out by getting a vector that is perpendicular to the given vector.
	V3D n(r, 0, 0);
	V3D dir = axis * direction.length();

	glBegin(GL_TRIANGLE_STRIP);

	//now, we we'll procede by rotating the vector n around v, and create the cylinder that way.
	for (int i = 0;i <= nPoints;i++) {
		V3D p = n.rotate(2 * i*PI / nPoints, axis);
		V3D normal = p.unit();
		glNormal3d(normal[0], normal[1], normal[2]);
		P3D p1 = org + p;
		glVertex3d(p1[0], p1[1], p1[2]);
		P3D p2 = org + dir + p;
		glVertex3d(p2[0], p2[1], p2[2]);
	}
	glEnd();

	//now we'll draw a half sphere...
	P3D p, q;
	double angle = PI / nPoints;
	p[0] = r*cos(-PI / 2);
	p[1] = r*sin(-PI / 2);
	p[2] = 0;
	for (int i = nPoints / 2;i >= 0;i--) {
		glBegin(GL_QUAD_STRIP);
		q[0] = r*cos(-i*angle);
		q[1] = r*sin(-i*angle);
		q[2] = 0;
		//make sure we compute the normal as well as the node coordinates
		n = V3D(p).toUnit();
		glNormal3d(n[0], n[1], n[2]);
		glVertex3d(org[0] + p[0], org[1] + p[1], org[2] + p[2]);
		//make sure we compute the normal as well as the node coordinates
		n = V3D(q).toUnit();
		glNormal3d(n[0], n[1], n[2]);
		glVertex3d(org[0] + q[0], org[1] + q[1], org[2] + q[2]);

		for (int j = 0;j <= nPoints;j++) {
			//make sure we compute the normal as well as the node coordinates
			V3D v = V3D(q[0] * cos(2 * j * angle), q[1], q[0] * sin(2 * j * angle));
			n = v.unit();
			glNormal3d(n[0], n[1], n[2]);
			glVertex3d(org[0] + v[0], org[1] + v[1], org[2] + v[2]);

			//make sure we compute the normal as well as the node coordinates
			v = V3D(p[0] * cos(2 * j * angle), p[1], p[0] * sin(2 * j * angle));
			n = v.unit();
			glNormal3d(n[0], n[1], n[2]);
			glVertex3d(org[0] + v[0], org[1] + v[1], org[2] + v[2]);
		}
		p = q;
		glEnd();
	}

	org = P3D() + dir;

	//	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	//and now draw the other half - horrible, quick research code...

	p[0] = r*cos(PI / 2);
	p[1] = r*sin(PI / 2);
	p[2] = 0;
	for (int i = nPoints / 2;i >= 0;i--) {
		glBegin(GL_QUAD_STRIP);
		q[0] = r*cos(i*angle);
		q[1] = r*sin(i*angle);
		q[2] = 0;
		//make sure we compute the normal as well as the node coordinates
		n = V3D(p).toUnit();
		glNormal3d(n[0], n[1], n[2]);
		glVertex3d(org[0] + p[0], org[1] + p[1], org[2] + p[2]);
		//make sure we compute the normal as well as the node coordinates
		n = V3D(q).toUnit();
		glNormal3d(n[0], n[1], n[2]);
		glVertex3d(org[0] + q[0], org[1] + q[1], org[2] + q[2]);

		for (int j = 0;j <= nPoints;j++) {
			//make sure we compute the normal as well as the node coordinates
			V3D v = V3D(q[0] * cos(2 * j * angle), q[1], q[0] * sin(2 * j * angle));
			n = v.unit();
			glNormal3d(n[0], n[1], n[2]);
			glVertex3d(org[0] + v[0], org[1] + v[1], org[2] + v[2]);

			//make sure we compute the normal as well as the node coordinates
			v = V3D(p[0] * cos(2 * j * angle), p[1], p[0] * sin(2 * j * angle));
			n = v.unit();
			glNormal3d(n[0], n[1], n[2]);
			glVertex3d(org[0] + v[0], org[1] + v[1], org[2] + v[2]);
		}
		p = q;
		glEnd();
	}

	glPopMatrix();
}

// draws a capsule from origin to origin + direction
void drawCapsule(const P3D& start, const P3D& end, double r, int nPoints) {
	drawCapsule(start, V3D(start, end), r, nPoints);
}

// draws a cylinder from origin to origin + direction
void drawCylinder(const P3D& origin, const V3D& direction, double r, int nPoints, bool cappedEnds) {
	if (IS_ZERO(direction.length()) || IS_ZERO(r)) return;

	//we'll start out by getting a vector that is perpendicular to the given vector.
	V3D n, t, axis = direction.unit();

	axis.getOrthogonalVectors(n, t);

	n = n.unit() * r;

	glBegin(GL_TRIANGLE_STRIP);
	//now, we we'll procede by rotating the vector n around v, and create the cylinder that way.
	for (int i = 0;i <= nPoints;i++) {
		V3D p = n.rotate(2 * i*PI / nPoints, -axis);
		V3D normal = p.unit();
		glNormal3d(normal[0], normal[1], normal[2]);
		P3D p1 = origin + p;
		glVertex3d(p1[0], p1[1], p1[2]);
		P3D p2 = origin + direction + p;
		glVertex3d(p2[0], p2[1], p2[2]);
	}
	glEnd();

/*
	glBegin(GL_LINES);
	//now, we we'll procede by rotating the vector n around v, and create the cylinder that way.
	for (int i = 0;i <= nPoints;i++) {
		V3D p = n.rotate(2 * i*PI / nPoints, axis);
		V3D normal = p.unit();
		P3D p1 = origin + p;
		P3D p2 = origin + direction + p;

		glVertex3d(p1[0], p1[1], p1[2]);
		glVertex3d(p1[0] + normal[0]*0.05, p1[1] + normal[1]*0.05, p1[2] + normal[2]*0.05);

		glVertex3d(p2[0], p2[1], p2[2]);
		glVertex3d(p2[0] + normal[0]*0.05, p2[1] + normal[1]*0.05, p2[2] + normal[2]*0.05);
	}
	glEnd();
*/
	if (cappedEnds) {
/*
		glBegin(GL_LINES);
		glVertex3d(origin[0], origin[1], origin[2]);
		glVertex3d(origin[0]-axis[0]*0.05, origin[1] - axis[1]*0.05, origin[2] - axis[2]*0.05);
		glEnd();
*/
		glBegin(GL_POLYGON);
		glNormal3d(-axis[0], -axis[1], -axis[2]);
		//now, we we'll procede by rotating the vector n around v, and create the cylinder that way.
		for (int i = 0;i < nPoints;i++) {
			P3D p1 = origin + n.rotate(2 * i*PI / nPoints, -axis);
			glVertex3d(p1[0], p1[1], p1[2]);
			P3D p2 = origin + n.rotate(2 * (i + 1)*PI / nPoints, -axis);
			glVertex3d(p2[0], p2[1], p2[2]);
		}
		glEnd();
/*
		glBegin(GL_LINES);
		glVertex3d(origin[0] + direction[0], origin[1]+ direction[1], origin[2] + direction[2]);
		glVertex3d(origin[0] + direction[0] + axis[0]*0.05, origin[1] + direction[1] + axis[1]*0.05, origin[2] + direction[2] + axis[2]*0.05);
		glEnd();
*/
		glBegin(GL_POLYGON);
		glNormal3d(axis[0], axis[1], axis[2]);
		//now, we we'll procede by rotating the vector n around v, and create the cylinder that way.
		for (int i = 0;i < nPoints;i++) {
			P3D p1 = origin + direction + n.rotate(2 * i*PI / nPoints, axis);
			glVertex3d(p1[0], p1[1], p1[2]);
			P3D p2 = origin + direction + n.rotate(2 * (i + 1)*PI / nPoints, axis);
			glVertex3d(p2[0], p2[1], p2[2]);
		}
		glEnd();
	}
}


/**
This method draws a cone of radius r, along the vector dir, with the center of its base at org.
*/
void drawCone(P3D org, V3D v, double r, int nrPoints) {
	//we'll start out by getting a vector that is perpendicular to the given vector.
	V3D n;
	V3D axis = v;
	axis.toUnit();
	int i;
	//try to get a vector that is not colinear to v.
	if (v[0] != 0 || v[1] != 0)
		n = V3D(v[0], v[1], v[2] + 1);
	else if (v[1] != 0 || v[2] != 0)
		n = V3D(v[0], v[1] + 1, v[2]);
	else
		n = V3D(v[0] + 1, v[1], v[2]);

	n = n.cross(v);

	if (IS_ZERO(v.length()) || IS_ZERO(n.length()))
		return;
	(n.toUnit()) *= r;
	glBegin(GL_TRIANGLE_FAN);

	P3D p2 = org + v;
	glNormal3d(axis[0], axis[1], axis[2]);
	glVertex3d(p2[0], p2[1], p2[2]);


	//now, we we'll procede by rotating the vector n around v, and creating the cone that way.
	for (i = 0;i <= nrPoints;i++) {
		V3D p = n.rotate(2 * i*PI / nrPoints, axis);
		V3D normal = p.unit();
		glNormal3d(normal[0], normal[1], normal[2]);
		P3D p1 = org + p;
		glVertex3d(p1[0], p1[1], p1[2]);
	}
	glEnd();


	//now we need to draw the bottom of the cone.
	glBegin(GL_POLYGON);

	//now, we we'll procede by rotating the vector n around v, and creating the cone that way.
	for (i = 0;i <= nrPoints;i++) {
		V3D p = n.rotate(2 * i*PI / nrPoints, axis);
		V3D normal = p.unit();
		glNormal3d(normal[0], normal[1], normal[2]);
		P3D p1 = org + p;
		glVertex3d(p1[0], p1[1], p1[2]);
	}
	glEnd();

}

// draws an arrow from start to end
void drawArrow(const P3D& start, const P3D& end, double r, int nPoints) {
	drawCylinder(start, start + V3D(start, end) * 0.75, r, nPoints);
	P3D arrowStart = start + V3D(start, end) * 0.75;
	P3D arrowEnd = start + V3D(start, end) * 1.0;
	drawCone(arrowStart, V3D(arrowStart, arrowEnd), 1.5*r);
}


// draws a cylinder from origin to origin + direction
void drawCylinder(const P3D& start, const P3D& end, double r, int nPoints, bool cappedEnds) {
	drawCylinder(start, V3D(start, end), r, nPoints, cappedEnds);
}

// draw an axis-aligned box defined by the two corner points
void drawBox(const P3D& min, const P3D& max) {
	//now draw the cube that is defined by these two points...
	glBegin(GL_QUADS);
	glNormal3d(0, 0, -1);
	glVertex3d(min[0], min[1], min[2]);
	glVertex3d(min[0], max[1], min[2]);
	glVertex3d(max[0], max[1], min[2]);
	glVertex3d(max[0], min[1], min[2]);

	glNormal3d(0, 0, 1);
	glVertex3d(min[0], min[1], max[2]);
	glVertex3d(max[0], min[1], max[2]);
	glVertex3d(max[0], max[1], max[2]);
	glVertex3d(min[0], max[1], max[2]);

	glNormal3d(0, -1, 0);
	glVertex3d(min[0], min[1], min[2]);
	glVertex3d(max[0], min[1], min[2]);
	glVertex3d(max[0], min[1], max[2]);
	glVertex3d(min[0], min[1], max[2]);

	glNormal3d(0, 1, 0);
	glVertex3d(min[0], max[1], min[2]);
	glVertex3d(min[0], max[1], max[2]);
	glVertex3d(max[0], max[1], max[2]);
	glVertex3d(max[0], max[1], min[2]);

	glNormal3d(-1, 0, 0);
	glVertex3d(min[0], min[1], min[2]);
	glVertex3d(min[0], min[1], max[2]);
	glVertex3d(min[0], max[1], max[2]);
	glVertex3d(min[0], max[1], min[2]);

	glNormal3d(1, 0, 0);
	glVertex3d(max[0], min[1], min[2]);
	glVertex3d(max[0], max[1], min[2]);
	glVertex3d(max[0], max[1], max[2]);
	glVertex3d(max[0], min[1], max[2]);

	glEnd();

}


void generateMenuEnumFromFileList(const std::string& menuName, const DynamicArray<std::string>& itemList) {
	std::string enumString;
	int index = 0;
	enumString += menuName;
	enumString += " enum='";
	for (auto it = itemList.begin(); it != itemList.end(); ++it) {
		if (it != itemList.begin())
			enumString += ',';
		enumString += std::to_string(index);
		enumString += " {";
		enumString += it->substr(it->find_last_of('\\') + 1);
		enumString += "}";
		index += 1;
	}
	enumString += '\'';

	//	Logger::consolePrint("--> %s\n", enumString.c_str());

	// Defining new enum type
	TwDefine(enumString.c_str());
}

void drawMOIApproximation(const Matrix3x3 MOI, double mass) {
	Eigen::EigenSolver<Matrix3x3> eigenvalueSolver(MOI);

	Eigen::Vector3cd principleMomentsOfInertia = eigenvalueSolver.eigenvalues();

	assert(IS_ZERO(principleMomentsOfInertia[0].imag()) && IS_ZERO(principleMomentsOfInertia[1].imag()) && IS_ZERO(principleMomentsOfInertia[1].imag()));

	Eigen::Matrix3cd V = eigenvalueSolver.eigenvectors();

	double Ixx = principleMomentsOfInertia[0].real(); // = m(y2 + z2)/12
	double Iyy = principleMomentsOfInertia[1].real(); // = m(z2 + x2)/12
	double Izz = principleMomentsOfInertia[2].real(); // = m(y2 + x2)/12

	double x = sqrt((Iyy + Izz - Ixx) * 6 / mass);
	double y = sqrt((Izz + Ixx - Iyy) * 6 / mass);
	double z = sqrt((Ixx + Iyy - Izz) * 6 / mass);

	P3D pmin(-x / 2, -y / 2, -z / 2), pmax(x / 2, y / 2, z / 2);

	if (V.determinant().real() < 0.0) {
		V(0, 2) *= -1;
		V(1, 2) *= -1;
		V(2, 2) *= -1;
	}
	assert(IS_ZERO(abs(V.determinant().real() - 1.0)) && "Rotation matrices have a determinant which is equal to 1.0!");

	double glmatrix[16]; memset(glmatrix, 0, sizeof(glmatrix));

	for (uint i = 0; i < 3; i++)
		for (uint j = 0; j < 3; j++)
		glmatrix[j * 4 + i] = V(i, j).real();
	glmatrix[15] = 1.0;

	glPushMatrix();
	glMultMatrixd(glmatrix);
	drawBox(pmin, pmax);
	glPopMatrix();
}

void checkOGLErrors() {
	int glErr = glGetError();
//	if (glErr != GL_NO_ERROR)
//		Logger::consolePrint("GL error: %d %s\n", glErr, gluErrorString(glErr));
}


