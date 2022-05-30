#include "GLUtils.h"
#include "RotateWidgetV2.h"
#include <MathLib/MathLib.h>
#include <Utils/Logger.h>
#include "GlobalMouseState.h"

RotateWidgetV2::RotateWidgetV2(uint axes) {
	if (axes & AXIS_X)
		rAxes.push_back(RotateWidgetV2Torus(this, V3D(1, 0, 0)));
	if (axes & AXIS_Y)
		rAxes.push_back(RotateWidgetV2Torus(this, V3D(0, 1, 0)));
	if (axes & AXIS_Z)
		rAxes.push_back(RotateWidgetV2Torus(this, V3D(0, 0, 1)));
}

RotateWidgetV2::~RotateWidgetV2() {
}

void RotateWidgetV2::draw() {
	if (visible) {
		glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_CURRENT_BIT);
		glEnable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);
		glPushMatrix();

		//to keep the size of the arrows roughly the same on screen, mess around with the scale...
		P3D camPos = getRayFromScreenCoords(0, 0).origin;
		double z = V3D(camPos, pos).length();
		scale = 50 * z / 1000;

		// Get and set the transformation matrix
		glTranslated(pos[0], pos[1], pos[2]);
		glColor3d(1, 1, 0);
		drawSphere(P3D(0, 0, 0), scale * 0.2);

		for (auto it = rAxes.begin(); it != rAxes.end(); ++it)
			it->draw();

		glPopMatrix();
		glPopAttrib();
		glColor3d(1, 1, 1);
	}
}

bool RotateWidgetV2::isPicked() {
	for (uint i = 0; i < rAxes.size(); i++) {
		if (rAxes[i].isPicked)
			return true;
	}
	return false;
}

//triggered when mouse moves
bool RotateWidgetV2::onMouseMoveEvent(double xPos, double yPos) {
	if (visible == false) return false;
	Ray mouseRay = getRayFromScreenCoords(xPos, yPos);

	if (GlobalMouseState::lButtonPressed) {
		for (auto it = rAxes.begin(); it != rAxes.end(); ++it)
			if (it->isPicked) {
				double newTorusAngle = it->getAngleOnTorusForRay(mouseRay);
	//			Logger::consolePrint("TWIST: %lf %lf\n", newTorusAngle, torus->lastPickedAngle);
				double angleOffset = newTorusAngle - it->lastPickedAngle;
				if (angleOffset > 2 * PI) angleOffset -= 2 * PI;
				if (angleOffset < -2 * PI) angleOffset += 2 * PI;
				//note: the twist angle is not bounded on either side, so it can increase arbitrarily...
				orientation = getRotationQuaternion(angleOffset, it->rotationAxis) * orientation;
				it->lastPickedAngle = newTorusAngle;
				return true;
			}
	}
	else {
		active = false;
		for (auto it = rAxes.begin(); it != rAxes.end(); ++it)
			it->isPicked = false;
		for (auto it = rAxes.begin(); it != rAxes.end(); ++it) {
			it->pickWith(mouseRay);
			if (it->isPicked) {
				active = true;
				it->lastPickedAngle = it->getAngleOnTorusForRay(mouseRay);
				break;
			}
		}
		if (active) return true;
	}
	return false;
}


double RotateWidgetV2Torus::getAngleOnTorusForRay(const Ray& ray) {
	Plane plane(pWidget->pos, rotationAxis);
	P3D pickedPoint;
	ray.getDistanceToPlane(plane, &pickedPoint);
	//get the rotation angle, relative to the zero axis, that corresponds to the picked point...
	return (zeroRotationAxis).angleWith(V3D(pWidget->pos, pickedPoint), rotationAxis);
}


void RotateWidgetV2Torus::pickWith(const Ray& ray) {
	Plane plane(pWidget->pos, rotationAxis);
	P3D pickedPoint;
	ray.getDistanceToPlane(plane, &pickedPoint);
	double dist = V3D(pickedPoint - pWidget->pos).length();

	if (abs(dist - pWidget->scale * 1.5) < pWidget->scale / 5) {
		isPicked = true;
	}
	else
		isPicked = false;

}

void RotateWidgetV2Torus::draw() {
	double alpha = isPicked ? 3 : 1;
	glColor3d(rotationAxis[0] / 3 * alpha, rotationAxis[1] / 3 * alpha, rotationAxis[2] / 3 * alpha);
	drawTorus(P3D(0, 0, 0), rotationAxis, pWidget->scale * 1.5 * (1 + (1 - 1) / 20.0), pWidget->scale * (1 + (10*alpha - 10) / 20.0) / 40);

//	P3D p1 = P3D(0, 0, 0) + (zeroRotationAxis).rotate(0, rotationAxis) * pWidget->scale * 1.5 * (1 + (alpha - 1) / 20.0);
//	drawCylinder(P3D(0,0,0), p1, pWidget->scale * (1 + (alpha - 1) / 20.0) / 20);
}
