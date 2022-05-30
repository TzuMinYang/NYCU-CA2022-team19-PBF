#include "InteractiveWidget.h"
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

InteractiveWidget::InteractiveWidget() {}

InteractiveWidget::~InteractiveWidget() {}


Ray InteractiveWidget::getRayFromScreenCoords(double x, double y) {
	Ray theRay;

	//get the modelview matrix
	double modelviewMat[16];
	glGetDoublev(GL_MODELVIEW_MATRIX, modelviewMat);

	//and get the projection matrix
	double projectionMat[16];
	glGetDoublev(GL_PROJECTION_MATRIX, projectionMat);

	//and the viewport info
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	if (viewport[2] <= 0 || viewport[3] <= 0)
		return Ray();

	P3D p1, p2;
	gluUnProject(x, viewport[3] - y - 1, 0, modelviewMat, projectionMat, viewport, &p1[0], &p1[1], &p1[2]);
	gluUnProject(x, viewport[3] - y - 1, 1, modelviewMat, projectionMat, viewport, &p2[0], &p2[1], &p2[2]);

	theRay.origin = p1;
	theRay.direction = V3D(p1, p2).unit();

	return theRay;
}
