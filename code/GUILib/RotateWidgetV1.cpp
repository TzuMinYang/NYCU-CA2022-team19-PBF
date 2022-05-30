#include "GLUtils.h"
#include "RotateWidgetV1.h"
#include <MathLib/MathLib.h>
#include <Utils/Logger.h>
#include "GlobalMouseState.h"

RotateWidgetV1::RotateWidgetV1() {
	axis = new RotateWidgetV1Axis(this);
	torus = new RotateWidgetV1Torus(this);
	pos = P3D(0, 0, 0);
}

RotateWidgetV1::~RotateWidgetV1() {
}

void RotateWidgetV1::draw() {
	if (visible) {
		glPushAttrib(GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_CURRENT_BIT);
		glEnable(GL_LIGHTING);
		glDisable(GL_TEXTURE_2D);

//		glDisable(GL_DEPTH_TEST);
		glPushMatrix();

		//to keep the size of the arrows roughly the same on screen, mess around with the scale...
		P3D camPos = getRayFromScreenCoords(0, 0).origin;
		double z = V3D(camPos, pos).length();
		scale = 50 * z / 1000;

		// Get and set the transformation matrix
		glTranslated(pos[0], pos[1], pos[2]);

		axis->draw();
		torus->draw();

		glPopMatrix();
		glPopAttrib();

		glColor3d(1, 1, 1);
	}
}

bool RotateWidgetV1::isPicked() {
	if (axis->isPicked || torus->isPicked)
		return true;
	else
		return false;
}

//triggered when mouse moves
bool RotateWidgetV1::onMouseMoveEvent(double xPos, double yPos) {
	if (visible == false) return false;
	V3D mouseMoveOffset(GlobalMouseState::mouseMoveX, GlobalMouseState::mouseMoveY);
	Ray mouseRay = getRayFromScreenCoords(xPos, yPos);

	if (GlobalMouseState::lButtonPressed) {
		P3D newPickedPoint; 
		double dist = mouseRay.getDistanceToPoint(pos, &newPickedPoint);
		double axisLen = scale * 1.5 * (1 + 1.0 / 20.0);

		if (axis->isPicked) {
			//V3D overallOffset = V3D(pickedPoint, newPickedPoint) / scale;
			//V3D widgetAxis = axis->getAxis() + overallOffset;
			V3D widgetAxis;
			if (dist > axisLen)
				widgetAxis = V3D(pos, newPickedPoint);
			else {
				//find the point that is on the ray, the unit sphere centered at pos, and closest to the current configuration...
				double offset = sqrt(axisLen*axisLen - dist*dist);
				P3D p1 = newPickedPoint + mouseRay.direction * offset;
				P3D p2 = newPickedPoint + mouseRay.direction * -offset;
				P3D pOriginal = pos + axis->getAxis();
//				Logger::consolePrint("dist: %lf, axisLen: %lf, dir: %lf, lenghts: %lf %lf\n", dist, axisLen, mouseRay.direction.length(), V3D(pos, p1).length(), V3D(pos, p2).length());

				double d1 = V3D(pOriginal, p1).length();
				double d2 = V3D(pOriginal, p2).length();
				if ( d1 < d2)
					widgetAxis = V3D(pos, p1);
				else
					widgetAxis = V3D(pos, p2);
			}

			widgetAxis.normalize();
			axis->setRotationToMatch(widgetAxis);
			return true;
		}

		if (torus->isPicked) {
			double newTorusAngle = torus->getAngleOnTorusForRay(mouseRay);
//			Logger::consolePrint("TWIST: %lf %lf\n", newTorusAngle, torus->lastPickedAngle);
			double angleOffset = newTorusAngle - torus->lastPickedAngle;
			if (angleOffset > 2 * PI) angleOffset -= 2 * PI;
			if (angleOffset < -2 * PI) angleOffset += 2 * PI;
			//note: the twist angle is not bounded on either side, so it can increase arbitrarily...
			torus->twistAngle += angleOffset;
			torus->lastPickedAngle = newTorusAngle;
			return true;
		}
	}
	else {

		active = false;
		axis->pickWith(mouseRay);
		if (axis->isPicked)
			active = true;

		torus->pickWith(mouseRay);
		if (torus->isPicked) {
			active = true;
			torus->lastPickedAngle = torus->getAngleOnTorusForRay(mouseRay);
		}

		if (active) return true;
	}
	return false;
}

void RotateWidgetV1Axis::pickWith(const Ray& ray) {
	P3D closestPtToRay;
	// Get the distance between the ray and the axis
	double dist = ray.getDistanceToSegment(pWidget->pos + getAxis() * pWidget->scale *0.5, pWidget->pos + getAxis() * pWidget->scale * 2, &closestPtToRay);

	if (dist < pWidget->scale / 5) {
		isPicked = true;
	}
	else
		isPicked = false;
}

void RotateWidgetV1Axis::draw() {
	double alpha = isPicked ? 2 : 1;
	glColor3d(1.0 / 2 * alpha, 1.0 / 2 * alpha, 1.0 / 2 * alpha);
	drawCylinder(P3D(), P3D() + getAxis()*pWidget->scale * 1.5 * (1 + (alpha - 1) / 20.0), pWidget->scale * (1 + (alpha - 1) / 20.0) / 20);
	drawSphere(P3D() + getAxis()*pWidget->scale * 1.5 * (1 + (alpha - 1) / 20.0), pWidget->scale * (1 + (alpha - 1) / 20.0) / 10);
}

double RotateWidgetV1Torus::getAngleOnTorusForRay(const Ray& ray) {
	Plane plane(pWidget->pos, pWidget->axis->getAxis());
	P3D pickedPoint;
	ray.getDistanceToPlane(plane, &pickedPoint);
	//get the rotation angle, relative to the zero axis, that corresponds to the picked point...
	return (pWidget->axis->rot * zeroAxis).angleWith(V3D(pWidget->pos, pickedPoint), pWidget->axis->getAxis());
}


void RotateWidgetV1Torus::pickWith(const Ray& ray) {
	Plane plane(pWidget->pos, pWidget->axis->getAxis());
	P3D pickedPoint;
	ray.getDistanceToPlane(plane, &pickedPoint);
	double dist = V3D(pickedPoint - pWidget->pos).length();

	if (abs(dist - pWidget->scale * 1.5) < pWidget->scale / 5) {
		isPicked = true;
	}
	else
		isPicked = false;
}

void RotateWidgetV1Torus::draw() {
	double alpha = isPicked ? 2 : 1;
	glColor3d(1.0 / 2 * alpha, 1.0 / 2 * alpha, 1.0 / 2 * alpha);
	drawTorus(P3D(0, 0, 0), pWidget->axis->getAxis(), pWidget->scale * 1.5 * (1 + (alpha - 1) / 20.0), pWidget->scale * (1 + (alpha - 1) / 20.0) / 20);

	P3D p1 = P3D(0, 0, 0) + (pWidget->axis->rot * zeroAxis).rotate(twistAngle, pWidget->axis->getAxis()) * pWidget->scale * 1.5 * (1 + (alpha - 1) / 20.0);
	drawCylinder(P3D(0,0,0), p1, pWidget->scale * (1 + (alpha - 1) / 20.0) / 20);
}
